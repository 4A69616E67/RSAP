import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

import com.github.SnowFlakes.File.CommonFile.CommonFile;
import com.github.SnowFlakes.System.Qsub;
import com.github.SnowFlakes.unit.Opts;
import com.github.SnowFlakes.unit.Parameter;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

public class Main {
    public File InputDir;
    public File OutputDir = new File("RSAP-Out");
    public String Prefix = "RSAP-Out";
    public String[] SampleName;
    public String[] GroupName;
    public File GenomeFile;
    public File IndexPrefix;
    public File AnnotationFile;
    public int thread = 4;
    public String CM = "HTSeq";
    public String attribute = "gene_name";
    public boolean PBS = false;

    public Main() {
    }

    public Main(File inputDir, File genomeFile) {
        InputDir = inputDir;
        GenomeFile = genomeFile;
    }

    public static void main(String[] args) throws Exception {
        Options Argument = new Options();
        Argument.addOption(Option.builder("d").longOpt("dir").hasArg().desc("input dir, dir format: {input dir}/{sample name}/{sample data})").required().build());
        Argument.addOption(Option.builder("g").longOpt("genome").hasArg().desc("genome file or genome index prefix").required().build());
        Argument.addOption(Option.builder("sample").hasArgs().desc("samples name or a file include sample name, separated with a space (program will use dir name in {input dir} as sample if user don't set this option)").build());
        Argument.addOption(Option.builder("group").hasArgs().desc("group name or file include group name").build());
        Argument.addOption(Option.builder("p").longOpt("prefix").hasArg().desc("out prefix").build());
        Argument.addOption(Option.builder("o").longOpt("out").hasArg().desc("out dir").build());
        Argument.addOption(Option.builder("t").longOpt("thread").hasArg().desc("thread for each sample (default 4)").build());
        Argument.addOption(Option.builder("anno").hasArg().argName("file").desc("annotation file (gtf file)").build());
        Argument.addOption(Option.builder("CM").longOpt("CountMethod").hasArg().argName("String").desc("\"HTSeq\" or \"featureCount\" or \"None\" (default \"HTSeq\")").build());
        Argument.addOption(Option.builder("attribute").hasArg().argName("String").desc("attribute in annotation file (\"gene_id\" or \"gene_name\" or other, default \"gene_name\")").build());
        Argument.addOption("pbs", "pbs", false, "commit command with pbs");
        if (args.length < 4) {
            new HelpFormatter().printHelp("java -jar " + Opts.JarFile.getName(), Argument, true);
            System.exit(1);
        }
        CommandLine cLine = new DefaultParser().parse(Argument, args);
        // ==============================================================================================
        Main main = new Main(Parameter.GetFileOpt(cLine, "d", null), Parameter.GetFileOpt(cLine, "g", null));
        main.SampleName = Parameter.GetStringOpts(cLine, "sample", null);
        main.GroupName = Parameter.GetStringOpts(cLine, "group", null);
        main.OutputDir = Parameter.GetFileOpt(cLine, "o", main.OutputDir);
        main.Prefix = Parameter.GetStringOpt(cLine, "p", main.Prefix);
        main.AnnotationFile = Parameter.GetFileOpt(cLine, "anno", null);
        main.thread = Parameter.GetIntOpt(cLine, "t", main.thread);
        main.CM = Parameter.GetStringOpt(cLine, "CM", main.CM);
        main.attribute = Parameter.GetStringOpt(cLine, "attribute", main.attribute);
        main.PBS = cLine.hasOption("pbs");
        main.run();

    }

    public boolean check() {
        boolean flag = true;
        if (!InputDir.isDirectory() || InputDir.listFiles().length == 0) {
            System.err.println("Uncorrect or empty input directory:\t" + InputDir);
            return false;
        }
        if (!GenomeFile.isFile()) {
            if (!GenomeFile.getParentFile().isDirectory()) {
                System.err.println("Uncorrect genome directory:\t" + GenomeFile.getParentFile());
                return false;
            }
            File[] list = GenomeFile.getParentFile().listFiles();
            int count = 0;
            for (File aFile : list) {
                if (aFile.getName().matches(GenomeFile.getName() + ".+")) {
                    count++;
                    IndexPrefix = GenomeFile;
                    GenomeFile = null;
                    break;
                }
            }
            if (count == 0) {
                System.err.println("Uncorrect genome index prefix:\t" + GenomeFile);
                return false;
            }

        }
        return flag;
    }

    public void run() throws IOException, InterruptedException {
        if (!check()) {
            return;
        }
        File[] SampleDirs = InputDir.listFiles();
        if (SampleName == null || SampleName.length == 0) {
            SampleName = new String[SampleDirs.length];
            for (int i = 0; i < SampleDirs.length; i++) {
                SampleName[i] = SampleDirs[i].getName();
            }
        }
        if (GroupName == null || GroupName.length == 0) {
            GroupName = new String[SampleName.length];
            for (int i = 0; i < GroupName.length; i++) {
                GroupName[i] = SampleName[i].replaceAll("_[^_]+$", "");
            }
        }
        System.out.println("Sample name:\t" + String.join(" ", SampleName));
        System.out.println("Group name:\t" + String.join(" ", GroupName));
        // -----------------------------------------create_index---------------------------------------------
        if (IndexPrefix == null) {
            String Command = "hisat2-build " + GenomeFile + " " + OutputDir + "/" + GenomeFile.getName();
            new com.github.SnowFlakes.System.CommandLine().run(Command, null, new PrintWriter(System.err));
            IndexPrefix = GenomeFile;
        }
        for (int i = 0; i < SampleDirs.length; i++) {
            File r1_fq = null, r2_fq = null;
            File[] list = SampleDirs[i].listFiles();
            for (int j = 0; j < list.length; j++) {
                if (list[j].getName().matches(".+_R1.+")) {
                    r1_fq = list[j];
                } else if (list[j].getName().matches(".+_R2.+")) {
                    r2_fq = list[j];
                }
            }
            if (r1_fq == null || r2_fq == null) {
                System.err.println("Error, can't find R1 or R2 file in " + SampleDirs[i]);
                System.exit(1);
            }
            File outDir = new File(OutputDir + "/" + SampleName[i]);
            outDir.mkdirs();
            CommonFile runFile = new CommonFile(outDir + "/run.sh");
            String commandLine = CreateCommand(r1_fq, r2_fq, IndexPrefix.toString(), outDir + "/" + Prefix + "_" + SampleName[i], thread, AnnotationFile);
            System.out.println("Alignment and filter: " + SampleName[i]);
            if (PBS) {
                Qsub qsub = new Qsub("1", thread, (long) 20e9, SampleName[i]);
                qsub.CreateSubmitFile(commandLine, runFile);
                qsub.run(runFile);
            } else {
                runFile.WriteOpen();
                runFile.WriteItem(commandLine);
                runFile.WriteClose();
                new com.github.SnowFlakes.System.CommandLine().run("sh " + runFile, new PrintWriter(System.out), new PrintWriter(System.err));
            }
        }

    }

    public String CreateCommand(File R1File, File R2File, String Index, String Prefix, int thread, File gtfFile) {
        thread = thread <= 0 ? 4 : thread;
        StringBuilder str = new StringBuilder();
        str.append("#alignment\n");
        File samFile = new File(Prefix + ".sam");
        str.append("hisat2 -p " + thread + " -x " + Index + " -1 " + R1File + " -2 " + R2File + " -S " + samFile + " 2> " + Prefix + ".align.log\n");
        str.append("#sam to bam and filter\n");
        File cleanBamFile = new File(Prefix + ".clean.bam");
        str.append("samtools view -@ " + (thread - 1) + " -Sbh -F 12 -o " + cleanBamFile + " " + samFile + "\n");
        str.append("#bam sort\n");
        File cleanSortBamFile = new File(Prefix + ".clean.sort.bam");
        str.append("samtools sort -@ " + (thread - 1) + " -o " + cleanSortBamFile + " " + cleanBamFile + "\n");
        str.append("#HTSeq count\n");
        if (gtfFile != null && gtfFile.exists()) {
            if (CM.equals("HTSeq")) {
                str.append("htseq-count -s no -r pos -i " + attribute + " -f bam " + cleanSortBamFile + " " + gtfFile + " 1> " + Prefix + "-HTSeq.count 2> " + Prefix + ".HTSeq.log\n");
            } else if (CM.equals("featureCount")) {
                str.append("featureCounts -T " + thread + " -p -g " + attribute + " -a " + gtfFile + " -o " + Prefix + "-featureCount.count " + cleanSortBamFile + "\n");
            }
        }
        return str.toString();
    }

}
