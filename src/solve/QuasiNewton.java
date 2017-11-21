package solve;


import func.QuadraForm;
import func.RealFunc;
import line.LineSearch;
import util.Matrix;
import util.Vector;

public class QuasiNewton extends Algorithm {
	private RealFunc f;
	private LineSearch s;
	private Matrix A0;
	private Matrix H0;
	private Vector x0;
	private double radius;
	@Override
	public void start(Vector x0) {
		this.A0 = Matrix.identity(x0.size());
		this.H0 = Matrix.identity(x0.size());
		this.x0 = x0;
		this.radius = 1.0;
		super.start(x0);
		
	}
	@Override
	protected void compute_next() throws EndOfIteration {
		// Creation du modele quadratique
		Vector OldIter = this.iter_vec;
		Vector OldGrad = this.f.grad(this.iter_vec);
		
		QuadraForm q = new QuadraForm(this.A0,OldGrad);
		
		
		/* Mise a jour des matricex*/
		Vector GradQuadra = q.grad(this.iter_vec);
		double alpha = this.s.search(this.iter_vec, GradQuadra);
		Vector PosTemp = this.iter_vec.sub(this.H0.mult(GradQuadra.leftmul(alpha)));
		if (q.eval(PosTemp)<q.eval(this.iter_vec))
		{
			double Rconf = this.iter_vec.sub(PosTemp).norm();
			if(Rconf <= this.radius)
			{
				this.iter_vec = PosTemp;
			}else {
				Vector GauchyDir = this.f.grad(this.iter_vec).leftmul(this.radius/this.f.grad(this.iter_vec).norm()).minus();
				alpha = this.s.search(this.iter_vec, GauchyDir);
				PosTemp = this.iter_vec.sub(this.H0.mult(GauchyDir.leftmul(alpha)));
				/* GAUCHY MOTHER FUCKER */
				/* Calcul intersection cercle de confiance et PosTemp de gauchy ? */
			}
			
		}
		
		
		
		Vector deltaG = this.f.grad(this.iter_vec).sub(OldGrad);
		Vector deltaX = this.iter_vec.sub(OldIter);
		/*Matrix CorrectiveMatrix*/
		double test = (deltaG.sub(this.A0.mult(deltaX))).scalar(deltaG.sub(this.A0.mult(deltaX)));
		double test2 = deltaX.scalar(deltaG.sub(this.A0.mult(deltaX)));
		
		this.A0 = this.A0.leftmul(test/test2);
		
		
		
	}

}
