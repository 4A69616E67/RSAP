import java.io.File;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

public class Main {
        private File InputDir;
        private File OutputDir = new File("RSAP-Out");
        private String[] SampleName;
        private String[] GroupName;
        private File GenomeFile;
        private String IndexPrefix = "RSAP-out";
        private File AnnotationFile;
        private int thread = 4;
        private String CM = "HTSeq";
        private String attribute = "gene_name";

        public static void main(String[] args) throws Exception {
                Options Argument = new Options();
                Argument.addOption(Option.builder("d").longOpt("dir").hasArg()
                                .desc("input dir, dir format: {input dir}/{sample name}/{sample data})").required()
                                .build());
                Argument.addOption(Option.builder("s").longOpt("sample").hasArgs().desc(
                                "samples name or a file include sample name, separated with a space (program will use dir name in {input dir} as sample if user don't set this option)")
                                .build());
                Argument.addOption(Option.builder("g").longOpt("group").hasArgs()
                                .desc("group name or file include group name").build());
                Argument.addOption(Option.builder("p").longOpt("prefix").hasArg().desc("out prefix").build());
                Argument.addOption(Option.builder("o").longOpt("out").hasArg().desc("out dir").build());
                Argument.addOption(Option.builder("t").longOpt("thread").hasArg().desc("thread (default 4)").build());
                Argument.addOption(Option.builder("anno").hasArg().argName("file").desc("annotation file (gtf file)")
                                .required().build());
                Argument.addOption(Option.builder("CM").longOpt("CountMethod").hasArg().argName("String")
                                .desc("\"HTSeq\" or \"featureCount\" or \"None\" (defaule \"HTSeq\")").build());
                Argument.addOption(Option.builder("attribute").hasArg().argName("String").desc(
                                "attribute in annotation file (\"gene_id\" or \"gene_name\" or other, defaule \"gene_name\")")
                                .build());
        }

        public String CreateCommand(File R1File, File R2File, String Index, String Prefix, int thread, File gtfFile) {
                thread = thread <= 0 ? 4 : thread;
                StringBuilder str = new StringBuilder();
                str.append("#alignment\n");
                str.append("hisat2 -p " + thread + " -x " + Index + " -1 " + R1File + " -2 " + R2File + " -S " + Prefix
                                + ".sam\n");
                str.append("#sam to bam and filter\n");
                str.append("samtools view @ " + (thread - 1) + " -Sbh -F 12 -o " + Prefix + ".clean.bam " + Prefix
                                + ".sam\n");
                str.append("#bam sort\n");
                str.append("samtools sort -@ " + (thread - 1) + " -o " + Prefix + ".clean.sort.bam " + Prefix
                                + ".clean.bam\n");
                str.append("#HTSeq count\n");
                if (CM.equals("HTSeq")) {
                        str.append("htseq-count -s no -r pos -i gene_id -f bam " + Prefix + ".clean.sort.bam " + gtfFile
                                        + " 1> " + Prefix + ".HTSeq.count 2> " + Prefix + ".HTSeq.log\n");
                } else if (CM.equals("featureCount")) {
                        str.append("featureCounts -T " + thread + " -p -g " + attribute + " -a " + gtfFile + " -o "
                                        + Prefix + "-featureCount.count " + Prefix + ".clean.sort.bam ");
                }
                return str.toString();
        }

}
