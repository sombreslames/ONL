package func;
import util.Matrix;
import util.Vector;

public class QuadraForm implements RealFunc {

	public Matrix Q;
	public Vector b;
	
	public QuadraForm(Matrix Q) {
		this.Q = new Matrix(Q);
		this.b = new Vector(Q.nb_cols());
		
	}
	public QuadraForm(Matrix Q,Vector b) {
		this.Q = new Matrix(Q);
		this.b = b;
	}
	@Override
	public double eval(Vector x) {
		
		/* TODO */
		/* le - b.scalar(x) est du tp2 */
		return x.scalar((this.Q.mult(x.leftmul(0.5)))) - this.b.scalar(x);
	}

	@Override
	public Vector grad(Vector x) {
		
		/* TODO */
		
		return this.Q.mult(x);
	}
	
	@Override
	public int dim() {
		return Q.nb_cols();
	}

}
