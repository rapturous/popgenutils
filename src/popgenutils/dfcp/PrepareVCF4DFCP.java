package popgenutils.dfcp;
/**
 * Classname
 * 
 * version 0.1
 *
 * Jul 6, 2015
 * 
 * The contents of this file are subject to the terms of the GNU
 * General Public License Version 3. For more information visit 
 * http://www.gnu.org/licenses/gpl.html
 */


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

/**
 * @author Ouroboros
 *
 */
public class PrepareVCF4DFCP {

	static String filename;
	static int window_size;
	static String command;
	static int overlap_size;
	static String output_dir;
	static double pmask=-1, pgeno=-1;
	static int amask=-1, ageno=-1, sgeno=-1;
	static int hetindex=-1;
	static String popmappingfile;
	static String popstokeep;
	
	String[] gt_parts;
	StringBuilder gt_string_builder;
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		CommandLineParser parser = new PosixParser();
		final Options options = new Options();

		options.addOption("help", false,"print this help.");
		
		// SPLIT COMMAND
		options.addOption("window", true,"used with [split] command; size of window");
		options.addOption("overlap", true,"used with [split] command; size of overlap in adjacent windows");
		
		// MASK COMMAND
		options.addOption("pmask", true,"used with [mask] command; percentage of alleles masked (made unknown)");
		options.addOption("amask", true,"used with [mask] command; absolute number of alleles masked (made unknown)");

		// MAKEGENO COMMAND
		options.addOption("pgeno", true,"used with [makegeno] command; percentage of individuals turned into genotypes");
		options.addOption("ageno", true,"used with [makegeno] command; absolute number of individuals turned into genotypes");
		options.addOption("sgeno", true,"used with [makegeno] command; specific individual turned into genotypes");

		// HET FILTER COMMAND
		options.addOption("hetindex", true,"used with [hetfilter] command; index of individual to filter the VCF's heterozygous sites");
		
		// POP FILTER COMMAND
		options.addOption("pops", true,"used with [popfilter] command; populations seperated by , to keep");
		options.addOption("mapfile", true,"used with [popfilter] command; sample to population map file");

		CommandLine cmdline = null;
		try {
			cmdline = parser.parse(options, args);
			args = cmdline.getArgs();
		} catch (ParseException e) {
			System.err.println(e.getMessage());
			HelpFormatter help = new HelpFormatter();
			help.printHelp("preparevcf.jar [options] [command] [filename] [output_dir]", options);
			System.exit(-1);
		}
		
		command = args[0];
		filename = args[1];
		output_dir = args[2] + "/";

		(new File(output_dir)).mkdirs();	
		
		
		if (cmdline.hasOption("window")) {
			window_size = Integer.parseInt(cmdline.getOptionValue("window"));
		}		
		if (cmdline.hasOption("overlap")) {
			overlap_size = Integer.parseInt(cmdline.getOptionValue("overlap"));
		}		
		if (cmdline.hasOption("pmask")) {
			pmask = Double.parseDouble(cmdline.getOptionValue("pmask"));
		}		
		if (cmdline.hasOption("amask")) {
			amask = Integer.parseInt(cmdline.getOptionValue("amask"));
		}		
		if (cmdline.hasOption("pgeno")) {
			pgeno = Double.parseDouble(cmdline.getOptionValue("pgeno"));
		}		
		if (cmdline.hasOption("ageno")) {
			ageno = Integer.parseInt(cmdline.getOptionValue("ageno"));
		}
		if (cmdline.hasOption("hetindex")) {
			hetindex = Integer.parseInt(cmdline.getOptionValue("hetindex"));
		}

		if (cmdline.hasOption("pops")) {
			popstokeep = cmdline.getOptionValue("pops");
		}
		if (cmdline.hasOption("mapfile")) {
			popmappingfile = cmdline.getOptionValue("mapfile");
		}
		
		
		
		PrepareVCF4DFCP preparer = new PrepareVCF4DFCP();
		
		switch (command.toLowerCase()) {
		case "split":
			preparer.split();
			break;
		case "mask":
			preparer.mask();
			break;
		case "makegeno":
			preparer.makegeno();
			break;

		case "hetfilter":
			preparer.filterhet();
			break;
			
		case "popfilter":
			preparer.filterpop();
			break;

		default:
			HelpFormatter help = new HelpFormatter();
			help.printHelp("preparevcf.jar [options] [command] [filename] [output_dir]", options);
			System.exit(-1);
			break;
		}
		
	}
	/**
	 * 
	 */
	private void filterpop() {
		Set<Integer> indices_to_keep = new HashSet<Integer>();
		Map<String,String> sample_to_pop = new HashMap<String,String>();
		Map<String,String> sample_to_superpop = new HashMap<String,String>();
		Set<String> pops_to_keep = new HashSet<String>();
		
		for (int i = 0; i < 9; i++) {
			indices_to_keep.add(i);
		}
		
		String[] popsparts = popstokeep.split(",");
		for (String pop : popsparts) {
			pops_to_keep.add(pop);
		}
		
		try (BufferedReader in = Files.newBufferedReader(Paths.get(popmappingfile),
				Charset.forName("UTF-8"))) {
			String line = null;
			while ((line = in.readLine()) != null) {
				String[] parts = line.split("\t");
				sample_to_pop.put(parts[0], parts[1]);
				sample_to_superpop.put(parts[0], parts[2]);
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		StringBuilder header = new StringBuilder();
		try (BufferedReader in = Files.newBufferedReader(Paths.get(filename),
				Charset.forName("UTF-8"))) {
			BufferedWriter out = null;

			String line = null;
			while ((line = in.readLine()) != null) {
			
				if(line.startsWith("#CHROM")) {
					//samples begin at 9
					out = Files.newBufferedWriter(Paths.get(output_dir + "/"+"popfilter_" + Paths.get(filename).getFileName()), Charset.forName("UTF-8"));
					out.write(header.toString());
					String[] parts = line.split("\t");
					for (int i = 9; i < parts.length; i++) {
						if(pops_to_keep.contains(sample_to_superpop.get(parts[i])))
								indices_to_keep.add(i);
					}
					out.write(parts[0]);
					for (int i = 1; i < parts.length; i++) {
						if(indices_to_keep.contains(i))
							out.write("\t"+parts[i]);
					}
					out.write(System.getProperty("line.separator"));
				} else if(line.startsWith("#")) {							
					header.append(line + System.getProperty("line.separator"));
				} else {
					// format at 8
					String[] parts = line.split("\t");
					out.write(parts[0]);
					for (int i = 1; i < parts.length; i++) {
						if(indices_to_keep.contains(i))
							out.write("\t"+parts[i]);
					}
					out.write(System.getProperty("line.separator"));

				}
			}
			out.close();
		} catch (IOException e) {
			System.err.println("could not read from file " + filename);
			e.printStackTrace();
		}
	}
	/**
	 * 
	 */
	private void filterhet() {
		StringBuilder header = new StringBuilder();
		try (BufferedReader in = Files.newBufferedReader(Paths.get(filename),
				Charset.forName("UTF-8"))) {
			BufferedWriter out = null;

			String line = null;
			while ((line = in.readLine()) != null) {
			
				if(line.startsWith("#CHROM")) {
					//samples begin at 9
					out = Files.newBufferedWriter(Paths.get(output_dir + "/"+hetindex+"_hetfilter" + Paths.get(filename).getFileName()), Charset.forName("UTF-8"));
					out.write(header.toString());
					out.write(line + System.getProperty("line.separator"));
				} else if(line.startsWith("#")) {							
					header.append(line + System.getProperty("line.separator"));
				} else {
					// format at 8
					String[] parts = line.split("\t");
					String[] parts2 = parts[8].split(":");
					int gt_pos = 0;
					for (int i = 1; i < parts2.length; i++) {
						if(parts2[i].equals("GT")) 
							gt_pos = i;
					}
					parts2 = parts[9+hetindex].split("\\||/");
					if(!parts2[0].equals(parts2[1]))
						out.write(line + System.getProperty("line.separator"));
				}
			}
			out.close();
		} catch (IOException e) {
			System.err.println("could not read from file " + filename);
			e.printStackTrace();
		}
	}
	/**
	 * 
	 */
	private void makegeno() {
		gt_string_builder = new StringBuilder();
		StringBuilder header = new StringBuilder();
		try (BufferedReader in = Files.newBufferedReader(Paths.get(filename),
				Charset.forName("UTF-8"))) {
			BufferedWriter out = null;
			BufferedWriter out_bgl_ref = null;
			BufferedWriter out_bgl_geno = null;

			String line = null;
			int number_of_individuals=0;
			Set<Integer> individuals_to_genotype = new HashSet<Integer>();
			while ((line = in.readLine()) != null) {
			
				if(line.startsWith("#CHROM")) {
					StringBuilder bglref_header = new StringBuilder();
					StringBuilder bglgeno_header = new StringBuilder();
					String[] parts = line.split("\t");
					//samples begin at 9
					number_of_individuals = parts.length-9;
					List<Integer> index_of_individuals = new ArrayList<Integer>();
					for (int i = 0; i < number_of_individuals; i++) {
						index_of_individuals.add(i);
					}
					Collections.shuffle(index_of_individuals);
					if(pgeno>=0) {
						// convert pgeno into ageno
						ageno=(int) Math.ceil(number_of_individuals*pgeno);
					}
					for (int i = 0; i < ageno; i++) {
						individuals_to_genotype.add(index_of_individuals.get(i));
					}

					bglref_header.append(parts[0]);
					bglgeno_header.append(parts[0]);
					for (int i = 1; i < 9; i++) {
						bglref_header.append("\t"+parts[i]);
						bglgeno_header.append("\t"+parts[i]);
					}
					for (int i = 9; i < parts.length; i++) {
						int index = i-9;
						if(individuals_to_genotype.contains(index)) {
							bglgeno_header.append("\t"+parts[i]);
						} else {
							bglref_header.append("\t"+parts[i]);
						}
					}
					out = Files.newBufferedWriter(Paths.get(output_dir + "/"+ageno+"_geno" + Paths.get(filename).getFileName()), Charset.forName("UTF-8"));
					out_bgl_ref = Files.newBufferedWriter(Paths.get(output_dir + "/"+ageno+"_geno_bglref" + Paths.get(filename).getFileName()), Charset.forName("UTF-8"));
					out_bgl_geno = Files.newBufferedWriter(Paths.get(output_dir + "/"+ageno+"_geno_bgl" + Paths.get(filename).getFileName()), Charset.forName("UTF-8"));
					out.write(header.toString());
					out_bgl_ref.write(header.toString());
					out_bgl_geno.write(header.toString());
					out.write(line + System.getProperty("line.separator"));
					out_bgl_geno.write(bglgeno_header + System.getProperty("line.separator"));
					out_bgl_ref.write(bglref_header + System.getProperty("line.separator"));
				} else if(line.startsWith("#")) {							
					header.append(line + System.getProperty("line.separator"));
				} else {
					// format at 8
					String[] parts = line.split("\t");
					String[] parts2 = parts[8].split(":");
					int gt_pos = 0;
					for (int i = 1; i < parts2.length; i++) {
						if(parts2[i].equals("GT")) 
							gt_pos = i;
					}
					out.write(parts[0]);
					out_bgl_geno.write(parts[0]);
					out_bgl_ref.write(parts[0]);
					for (int i = 1; i < parts.length; i++) {
						if(i>8) {
							parts2 = parts[i].split(":");
							if(gt_pos==0 && individuals_to_genotype.contains(i-9)) {
								out.write("\t"+convertToGT(parts2[0]));
								out_bgl_geno.write("\t"+convertToGT(parts2[0]));
							} else {
								out.write("\t"+parts2[0]);
								out_bgl_ref.write("\t"+parts2[0]);
							}
							for (int j = 1; j < parts2.length; j++) {
								if(gt_pos==j && individuals_to_genotype.contains(i-9)) {
									out.write(":"+convertToGT(parts2[i]));
									out_bgl_geno.write(":"+convertToGT(parts2[i]));
								} else {
									out.write(":"+parts2[i]);
									out_bgl_ref.write(":"+parts2[i]);
								}
							}
						} else {
							out.write("\t"+parts[i]);	
							out_bgl_ref.write("\t"+parts[i]);	
							out_bgl_geno.write("\t"+parts[i]);							
						}
					}
					out.write(System.getProperty("line.separator"));
					out_bgl_ref.write(System.getProperty("line.separator"));
					out_bgl_geno.write(System.getProperty("line.separator"));
				}
			}
			out.close();
			out_bgl_ref.close();
			out_bgl_geno.close();
		} catch (IOException e) {
			System.err.println("could not read from file " + filename);
			e.printStackTrace();
		}
	}
	/**
	 * @param string
	 * @return
	 */
	private String convertToGT(String gt_field) {
		gt_string_builder.setLength(0);
		gt_parts = gt_field.split("\\||/");
		gt_string_builder.append(gt_parts[0]);
		for (int i = 1; i < gt_parts.length; i++) {
			gt_string_builder.append("/" + gt_parts[i]);
		}
		return gt_string_builder.toString();
	}
	/**
	 * 
	 */
	private void mask() {
		StringBuilder header = new StringBuilder();
		try (BufferedReader in = Files.newBufferedReader(Paths.get(filename),
				Charset.forName("UTF-8"))) {
			BufferedWriter out = null;

			String line = null;
			int number_of_individuals=0;
			Set<Integer> individuals_to_genotype = new HashSet<Integer>();
			while ((line = in.readLine()) != null) {
			
				if(line.startsWith("#CHROM")) {
					//samples begin at 9
					number_of_individuals = line.split("\t").length-9;
					List<Integer> index_of_individuals = new ArrayList<Integer>();
					for (int i = 0; i < number_of_individuals; i++) {
						index_of_individuals.add(i);
					}
					Collections.shuffle(index_of_individuals);
					if(pmask>=0) {
						// convert pgeno into ageno
						amask=(int) Math.ceil(number_of_individuals*pmask);
					}
					for (int i = 0; i < amask; i++) {
						individuals_to_genotype.add(index_of_individuals.get(i));
					}
					out = Files.newBufferedWriter(Paths.get(output_dir + "/"+amask+"_mask" + Paths.get(filename).getFileName()), Charset.forName("UTF-8"));
					out.write(header.toString());
					out.write(line + System.getProperty("line.separator"));
				} else if(line.startsWith("#")) {							
					header.append(line + System.getProperty("line.separator"));
				} else {
					// format at 8
					String[] parts = line.split("\t");
					String[] parts2 = parts[8].split(":");
					int gt_pos = 0;
					for (int i = 1; i < parts2.length; i++) {
						if(parts2[i].equals("GT")) 
							gt_pos = i;
					}
					out.write(parts[0]);
					for (int i = 1; i < parts.length; i++) {
						if(i>8) {
							parts2 = parts[i].split(":");
							if(gt_pos==0  && individuals_to_genotype.contains(i-8)) {
								out.write("\t"+maskAlleles(parts2[0]));
							} else out.write(parts2[0]);
							for (int j = 1; j < parts2.length; j++) {
								if(gt_pos==j  && individuals_to_genotype.contains(i-8)) {
									out.write(":"+maskAlleles(parts2[i]));
								} else out.write(":"+parts2[i]);
							}
						} else {
							out.write("\t"+parts[i]);							
						}
					}
					out.write(System.getProperty("line.separator"));
				}
			}
		} catch (IOException e) {
			System.err.println("could not read from file " + filename);
			e.printStackTrace();
		}
	}
	/**
	 * @param string
	 * @return
	 */
	private String maskAlleles(String gt_field) {
		gt_string_builder.setLength(0);
		gt_parts = gt_field.split("\\||/");
		gt_string_builder.append(".");
		for (int i = 1; i < gt_parts.length; i++) {
			gt_string_builder.append("/.");
		}
		return gt_string_builder.toString();
	}
	/**
	 * 
	 */
	private void split() {
		StringBuilder header = new StringBuilder();
		LineBuilder lines = new LineBuilder(window_size);
		try (BufferedReader in = Files.newBufferedReader(Paths.get(filename),
				Charset.forName("UTF-8"))) {
			int cnt = 0;
			int filecnt = 0;
			String line = null;
			while ((line = in.readLine()) != null) {
				if(line.startsWith("#")) {
					header.append(line + System.getProperty("line.separator"));
				} else {
					cnt++;
					lines.addLine(line);
					if(cnt==window_size) {
						cnt=overlap_size;
						lines.printFile(header.toString(), output_dir+filename.substring(0, filename.indexOf('.'))+
								"_"+window_size+"_"+overlap_size+"_"+(filecnt++)+".vcf");
					}
				}
			}
		} catch (IOException e) {
			System.err.println("could not read from file " + filename);
			e.printStackTrace();
		}
	}
	
	private class LineBuilder {
		String[] lines;
		int startPos = 0;

		public LineBuilder(int line_count) {
			super();
			lines = new String[line_count];
		}
		
		public void addLine(String line) {
			lines[startPos]=line;
			startPos++;
			startPos=startPos%window_size;
		}
		
		public void printFile(String header, String filename) {
			
			try (BufferedWriter out = Files.newBufferedWriter(
					Paths.get(filename), Charset.forName("UTF-8"))) {
				out.write(header);

				for (int i = startPos; i < startPos+lines.length; i++) {
					out.write(lines[i%lines.length] + System.getProperty("line.separator"));
				}
			} catch (IOException e) {
				System.err.println("could not write to file " + filename);
				e.printStackTrace();
			}
		}
	}

}
