package solve;
import line.LineSearch;
import func.RealFunc;
import util.Vector;

/**
 * Basic steepest descent algorithm for unconstrained minimization problem.
 * 
 * @author Gilles Chabert
 */

public class SteepestDescent extends Algorithm {
	
	private RealFunc f;
	private LineSearch s;
	
	/**
	 * Build the algorithm
	 * 
	 * @param f function to minimize
	 * @param s underlying line search algorithm
	 */
	public SteepestDescent(RealFunc f, LineSearch s) {
		this.f = f;
		this.s = s;
	}

	
	/**
	 * Calculate the next iterate.
	 * 
	 * Precondition:  the current iterate is stored in "iter_vec".
	 * 
	 * Postcondition: "iter_vec" is a reference to the next iterate. It can 
	 *                be obtained by calling current_vector(). If there is no 
	 *                more iterate, an exception "EndOfIteration" is thrown.
	 */
	public void compute_next() throws EndOfIteration {
		/* Direction */
		Vector Direction = new Vector(this.iter_vec.size());
		Direction = this.f.grad(this.iter_vec).minus();
		/* Calcul de alpha avec la methode la secante*/
		
		Vector ak = new Vector(this.iter_vec.size());		
		ak = this.iter_vec;
		
		/* Prochaine iteration */
		Vector akNext = Direction.leftmul(this.s.search(ak, Direction));
		akNext = ak.add(akNext);

		/* TODO */
		
		iter_vec= akNext ;
		// TODO: set iter_vec to xk+1
	}
	
}
