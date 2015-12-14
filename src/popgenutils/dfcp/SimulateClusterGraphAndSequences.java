package popgenutils.dfcp;
/**
 * Classname
 * 
 * version 0.1
 *
 * Sep 30, 2015
 * 
 * The contents of this file are subject to the terms of the GNU
 * General Public License Version 3. For more information visit 
 * http://www.gnu.org/licenses/gpl.html
 */


import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * @author Ouroboros
 *
 */
public class SimulateClusterGraphAndSequences {

	static Random r;

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int min_size_blocks = Integer.parseInt(args[0]);
		int max_size_blocks = Integer.parseInt(args[1]);
		double alpha = Double.parseDouble(args[2]);
		int min_splits_between_blocks = Integer.parseInt(args[3]);
		int max_splits_between_blocks = Integer.parseInt(args[4]);
		int number_of_individuals = Integer.parseInt(args[5]);
		r = new Random();
		r.setSeed(Long.parseLong(args[6]));
	}

	private List<Set<Integer>> sampleCRP(double alpha, int number_of_individuals) {
		List<Set<Integer>> partitioning = new ArrayList<Set<Integer>>();

		for (int i = 0; i < number_of_individuals; i++) {
			if(partitioning.isEmpty()) {
				partitioning.add(new HashSet<Integer>());
				partitioning.get(0).add(i);
			} else {
				double random_d = r.nextDouble();
				for (Set<Integer> table : partitioning) {
					
				}
			}
		}

		return partitioning;
	}
}
