package solve;
import line.LineSearch;
import func.RealFunc;
import util.Vector;

/**
 * Non-linear variant of the conjugate gradients.
 * 
 * @author Gilles Chabert
 *
 */
public class ConjugateGradients extends Algorithm {

	/**
	 * Build the algorithm for a given function and
	 * with an underlying line search technique.
	 * 
	 * @param f the function
	 * @param s the line search algorithm
	 */
	/**
	 * First direction of the current search (given by "start")
	 */
	protected Vector d;
	/**
	 * The function to minimize.
	 */
	protected RealFunc f;
	/**
	 * Point of the intipos search (given by "start")
	 */
	protected Vector x0;
	protected LineSearch s;
	
	public ConjugateGradients(RealFunc f, LineSearch s) {
		this.f = f;
		this.s = s;
		
	
	}
	
	/**
	 * Start the iteration
	 */
	public void start(Vector x0) {
		this.x0 = x0;
		this.d  = this.f.grad(x0).minus();
		super.start(x0);

	}
	
	/**
	 * Calculate the next iterate.
	 * 
	 * (update iter_vec).
	 */
	public void compute_next() throws EndOfIteration {
		
		Vector gk = this.f.grad(iter_vec);
		/* calcul de la nouvelle pos*/
		double alpha = this.s.search(iter_vec, d);
		Vector PosTemp = iter_vec.add(this.d.leftmul(alpha));
		
		
		Vector gk1 = this.f.grad(PosTemp).minus();
		/* Calcul de la nouvelle direction */
		double Beta = Math.pow(gk1.norm(),2)/Math.pow(gk.norm(),2);
		this.d = gk1.add(this.d.leftmul(Beta));
		iter_vec= PosTemp;
	}
}
