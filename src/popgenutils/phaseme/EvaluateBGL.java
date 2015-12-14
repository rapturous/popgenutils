/**
 * Classname
 * 
 * version 0.1
 *
 * Aug 31, 2015
 * 
 * The contents of this file are subject to the terms of the GNU
 * General Public License Version 3. For more information visit 
 * http://www.gnu.org/licenses/gpl.html
 */
package popgenutils.phaseme;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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

/**
 * @author Ouroboros
 *
 */
public class EvaluateBGL {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String bgl_output = args[0];
		String overall_vcf = args[1];


		Map<String,StringBuilder> individual_to_haplotype1 = new HashMap<String,StringBuilder>();
		Map<String,StringBuilder> individual_to_haplotype2 = new HashMap<String,StringBuilder>();
		Map<Integer,String> individual_map = new HashMap<Integer,String>();
		try (BufferedReader in = Files.newBufferedReader(Paths.get(bgl_output),
				Charset.forName("UTF-8"))) {
			StringBuilder bgl_results = new StringBuilder();

			String line = null;
			while ((line = in.readLine()) != null) {

				if(line.startsWith("#CHROM")) {
					String[] parts = line.split("\t");
					//samples begin at 9
					for (int i = 9; i < parts.length; i++) {
						individual_to_haplotype1.put(parts[i], new StringBuilder());
						individual_to_haplotype2.put(parts[i], new StringBuilder());
						individual_map.put(i,parts[i]);
					}
				} else if(line.startsWith("#")) {		

				} else {
					// format at 8
					String[] parts = line.split("\t");
					String[] parts2 = parts[8].split(":");
					int gt_pos = 0;
					for (int i = 1; i < parts2.length; i++) {
						if(parts2[i].equals("GT")) 
							gt_pos = i;
					}
					for (int i = 1; i < parts.length; i++) {
						if(i>8) {
							parts2 = parts[i].split(":");
							String[] GT = parts2[gt_pos].split("\\|");
							individual_to_haplotype1.get(individual_map.get(i)).append(GT[0]);
							individual_to_haplotype2.get(individual_map.get(i)).append(GT[1]);
						}
					}
				}
			}
		} catch (IOException e) {
			System.err.println("could not read from file " + bgl_output);
			e.printStackTrace();
		}

		Map<String,StringBuilder> individual_to_truehaplotype1 = new HashMap<String,StringBuilder>();
		Map<String,StringBuilder> individual_to_truehaplotype2 = new HashMap<String,StringBuilder>();
		try (BufferedReader in = Files.newBufferedReader(Paths.get(overall_vcf),
				Charset.forName("UTF-8"))) {

			individual_map.clear();
			String line = null;
			while ((line = in.readLine()) != null) {

				if(line.startsWith("#CHROM")) {
					String[] parts = line.split("\t");
					//samples begin at 9
					for (int i = 9; i < parts.length; i++) {
						individual_to_truehaplotype1.put(parts[i], new StringBuilder());
						individual_to_truehaplotype2.put(parts[i], new StringBuilder());
						individual_map.put(i,parts[i]);
					}
				} else if(line.startsWith("#")) {		

				} else {
					// format at 8
					String[] parts = line.split("\t");
					String[] parts2 = parts[8].split(":");
					int gt_pos = 0;
					for (int i = 1; i < parts2.length; i++) {
						if(parts2[i].equals("GT")) 
							gt_pos = i;
					}
					for (int i = 9; i < parts.length; i++) {
						if(i>8) {
							parts2 = parts[i].split(":");
							String[] GT = parts2[gt_pos].split("\\||/");
							individual_to_truehaplotype1.get(individual_map.get(i)).append(GT[0]);
							individual_to_truehaplotype2.get(individual_map.get(i)).append(GT[1]);
						}
					}
				}
			}
		} catch (IOException e) {
			System.err.println("could not read from file " + bgl_output);
			e.printStackTrace();
		}
		for (String ind : individual_to_haplotype1.keySet()) {
			System.err.println(ind);
			System.err.println(individual_to_haplotype1.get(ind).toString());
			System.err.println(individual_to_haplotype2.get(ind).toString());
			System.err.println(individual_to_truehaplotype1.get(ind).toString());
			System.err.println(individual_to_truehaplotype2.get(ind).toString());
			computeSwitchError(individual_to_haplotype1.get(ind).toString(),
					individual_to_haplotype2.get(ind).toString(),
					individual_to_truehaplotype1.get(ind).toString(),
					individual_to_truehaplotype2.get(ind).toString());
		}

	}


	private static double computeSwitchError(String hap1, String hap2, String hap1sol, String hap2sol) {
		int switches = 0;
		int non_matching_sites = 0;
		int hap1index = 0;
		int numHet = 0;
		for (int i = 0; i < hap1.length(); i++) {
			if(hap1sol.charAt(i)!=hap2sol.charAt(i)) {
				if(hap1sol.charAt(0)==hap1.charAt(i)) {
					hap1index=0;
					break;
				}
				if(hap2sol.charAt(0)==hap1.charAt(i)) {
					hap1index=1;
					break;
				}
			}
		}

		for (int i = 0; i < hap1.length(); i++) {
			byte[] genotype = null;
			genotype = new byte[] {Byte.valueOf(hap1sol.substring(i,i+1)),Byte.valueOf(hap2sol.substring(i,i+1))};
			if(genotype[0]!=genotype[1])
				numHet++;
			if(hap1index==0) {
				if(Byte.toString(genotype[0]).charAt(0)==hap1.charAt(i) &&
						Byte.toString(genotype[1]).charAt(0)==hap2.charAt(i)) {

				} else if(Byte.toString(genotype[1]).charAt(0)==hap1.charAt(i) &&
						Byte.toString(genotype[0]).charAt(0)==hap2.charAt(i)) {
					// switch
					switches++;
					hap1index=1;
				} else {
					non_matching_sites++;
				}
			} else {
				if(Byte.toString(genotype[1]).charAt(0)==hap1.charAt(i) &&
						Byte.toString(genotype[0]).charAt(0)==hap2.charAt(i)) {

				} else if(Byte.toString(genotype[0]).charAt(0)==hap1.charAt(i) &&
						Byte.toString(genotype[1]).charAt(0)==hap2.charAt(i)) {
					// switch
					switches++;
					hap1index=0;
				} else {
					non_matching_sites++;
				}
			}

		}
		System.err.println("switch errors:  " + switches);
		System.err.println("non matches: " + non_matching_sites);
		System.err.println("num sites: " + hap1.length());
		System.err.println("num het: " + numHet);
		return switches+non_matching_sites;
	}

}
