/**
 * Classname
 * 
 * version 0.1
 *
 * Oct 21, 2014
 * 
 * The contents of this file are subject to the terms of the GNU
 * General Public License Version 3. For more information visit 
 * http://www.gnu.org/licenses/gpl.html
 */
package popgenutils.phaseme;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * @author Ouroboros
 *
 */
public class EvalOutput {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		String filename = args[0];
		int num_runs = Integer.parseInt(args[1]);
		BufferedReader br = new BufferedReader(new FileReader(filename));
		int total1 = 0;
		int total2 = 0;
		int total3 = 0;
		for (int i = 0; i < num_runs; i++) {
			br.readLine();
			String[] real = new String[2];
			real[0] = br.readLine();
			real[1] = br.readLine();		
			// conditional phaser
			br.readLine();
			String[] computed = new String[2];
			computed[0] = br.readLine();
			computed[1] = br.readLine();	
			int dis = switchDistance(real,computed, false);
			int dis2 = switchDistance(real,computed, true);
			total1+=Math.min(dis, dis2);
			// joint phaser1
			br.readLine();
			computed[0] = br.readLine();
			computed[1] = br.readLine();	
			dis = switchDistance(real,computed, false);
			dis2 = switchDistance(real,computed, true);
			total2+=Math.min(dis, dis2);
			// joint phaser2
			br.readLine();
			computed[0] = br.readLine();
			computed[1] = br.readLine();	
			dis = switchDistance(real,computed, false);
			dis2 = switchDistance(real,computed, true);
			total3+=Math.min(dis, dis2);
		}

		System.out.println(((double)total1/(double)num_runs));
		System.out.println(((double)total2/(double)num_runs));
		System.out.println(((double)total3/(double)num_runs));
		
		br.close();
	}

	/**
	 * @param real
	 * @param computed
	 * @return
	 */
	private static int switchDistance(String[] real, String[] computed, boolean cis) {
		int dis = 0;
		for (int i = 0; i < real[0].length(); i++) {
			if(real[0].charAt(i)!=real[1].charAt(i)) {
				// het position
				if(cis) {
					if(real[0].charAt(i)!=computed[0].charAt(i)) {
						dis++;
						cis=false;
					}
				} else {
					if(real[0].charAt(i)!=computed[1].charAt(i)) {
						dis++;
						cis=true;
					}
				}
			}
		}	
		return dis;
	}


}
