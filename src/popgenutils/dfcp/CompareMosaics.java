package popgenutils.dfcp;


import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class CompareMosaics {

	public static void main(String[] args) throws IOException {

		if(args.length!=5) {
			System.err.println("args: number of haps, number of SNPs, mosaicfile1, mosaicfile2, number_compare");
			System.exit(-1);
		}
		int no_haps = Integer.parseInt(args[0]);
		int no_snps = Integer.parseInt(args[1]);
		String f1 = args[2];
		String f2 = args[3];
		int total_num_snps_to_sample = Integer.parseInt(args[4]);
		Set<Integer> randomSnps = new HashSet<Integer>();

		List<Integer> snp_nos = new ArrayList<Integer>();
		for (int i = 0; i < no_snps; i++) {
			snp_nos.add(i);
		}
		Collections.shuffle(snp_nos);
		for (int i = 0; i < total_num_snps_to_sample; i++) {
			randomSnps.add(snp_nos.get(i));
		}

		BufferedReader br = new BufferedReader(new FileReader(f1));
		BufferedReader br2 = new BufferedReader(new FileReader(f2));
		String line = null, line2 = null;
		String[] parts = null, parts2 = null;

		int[][] mosaic1 = new int[total_num_snps_to_sample][no_haps];
		int[][] mosaic2 = new int[total_num_snps_to_sample][no_haps];

		int lineno = 0;
		int snpno = 0;
		while((line=br.readLine())!=null) {
			line2 = br2.readLine();
			if(randomSnps.contains(lineno)) {
				parts = line.split(" ");
				parts2 = line2.split(" ");
				// labels don't matter, relationships among haplotypes do
				for (int i = 0; i < parts.length; i++) {
					mosaic1[snpno][i]=Integer.parseInt(parts[i]);
					mosaic2[snpno][i]=Integer.parseInt(parts2[i]);
				}
				snpno++;
			}
			lineno++;
		}
		double count_preserved = 0;
		double count_not_preserved = 0;
		double count_preserved2 = 0;
		double count_not_preserved2 = 0;
		for (int k = 0; k < total_num_snps_to_sample; k++) {
			for (int i = 0; i < no_haps; i++) {
				for (int j = i+1; j < no_haps; j++) {
					//Map<Integer,Set<Integer>> hap_to_hapsInGroup = new HashMap<Integer,Set<Integer>>();
					//Map<Integer,Set<Integer>> grpNumber_to_HapsInGroup = new HashMap<Integer,Set<Integer>>();					
					if(mosaic1[k][i]==mosaic1[k][j]) {
						if(mosaic2[k][i]==mosaic2[k][j]) {
							count_preserved++;
						} else count_not_preserved++;
					}
					if(mosaic2[k][i]==mosaic2[k][j]) {
						if(mosaic1[k][i]==mosaic1[k][j]) {
							count_preserved2++;
						} else count_not_preserved2++;
					}
				}
			}			
		}
		System.err.println(count_preserved);
		System.err.println(count_not_preserved);
		System.err.println(count_preserved2);
		System.err.println(count_not_preserved2);
		System.err.println(count_preserved/(count_preserved+count_not_preserved));
		System.err.println(count_preserved2/(count_preserved2+count_not_preserved2));
	}

}
