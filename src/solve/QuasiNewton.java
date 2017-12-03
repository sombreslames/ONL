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
	public QuasiNewton(RealFunc f,LineSearch s)
	{
		this.f=f;
		this.s=s;
	}
	@Override
	public void start(Vector x0) {
		this.A0 = Matrix.identity(x0.size());
		this.H0 = Matrix.identity(x0.size());
		this.x0 = x0;
		this.radius = DELTA_INIT;
		super.start(x0);
		
	}
	@Override
	protected void compute_next() throws EndOfIteration {
		Vector Direction = new Vector(this.iter_vec.size());
		Vector PosTemp 	= new Vector(this.iter_vec.size());
		double alpha   	= 0.0;
		double Adequacy = 0.0;
		double Rconf	= 0.0;
		Vector OldIter 	= this.iter_vec;
		Vector OldGrad 	= this.f.grad(this.iter_vec);
		
		// Creation du modele quadratique
		QuadraForm q = new QuadraForm(this.A0,OldGrad);
		
		Boolean GoodDirection = false;
		Boolean AlrightThen = false;//Booleen pour l'adequation
		Boolean GoGauchy = false;
		
		do {
			GoGauchy 		= false;
			GoodDirection 	= false;
			
			Direction = this.H0.mult(q.grad(this.iter_vec)).minus();
			alpha = this.s.search(this.iter_vec, Direction);
			PosTemp = this.iter_vec.add(Direction.leftmul(alpha));
			
			if ( f.eval(PosTemp) <= f.eval(this.iter_vec) )
			{
				Rconf = this.iter_vec.sub(PosTemp).norm();				
				if(Rconf > this.radius)
				{	
					GoGauchy = true;
					System.out.println("Descent direction not in the trust region\nRconf :"+Math.round(Rconf)+" radius :"+this.radius);
				}else {
					GoodDirection = true;
					System.out.println("Descent direction in the trust region");
				}
				
			}else {
				
				System.out.println("Not a descent direction");
				GoGauchy = true;
			}
			
			if(GoGauchy){
				System.out.println("So we compute gauchy : ");
				Vector GauchyDir = q.grad(this.iter_vec).leftmul(this.radius/q.grad(this.iter_vec).norm()).minus();
				alpha = this.s.search(this.iter_vec, GauchyDir);
				PosTemp = this.iter_vec.add(GauchyDir.leftmul(alpha));
				GoodDirection = true;
				/* Calcul intersection cercle de confiance et PosTemp de gauchy ? */
			}
			
			
			Adequacy = Math.abs((this.f.eval(OldIter)-this.f.eval(PosTemp))/(q.eval(OldIter)-q.eval(PosTemp)));
			AlrightThen = false;
			if(Adequacy >= GOOD_ADEQUACY)
			{
				this.radius *= DELTA_RATIO;
				this.iter_vec= PosTemp;
				System.out.println("Good adequacy : "+Adequacy);
				System.out.println("It is working !");
			}else if (Adequacy < POOR_ADEQUACY)
			{
				System.out.println("Bad adequacy : "+Adequacy);
				this.radius /= DELTA_RATIO;
				AlrightThen = true;
			}
		}while(AlrightThen && GoodDirection);
		System.out.println("Iteration "+this.current_iteration());
		System.out.println("f(k) : "+f.eval(OldIter)+"\nf(k+1) : "+f.eval(this.iter_vec));
		
		
		
		
		/* Mise a jour des matricex SR1*/
		/* DELTA COMPUTATION*/
		Vector deltaX 		= this.iter_vec.sub(OldIter);
		Vector deltaG 		= q.grad(this.iter_vec).sub(OldGrad);
		/* Error  */
		Vector ErrorKA		= deltaG.sub(this.A0.mult(deltaX));
		Vector ErrorKH		= deltaX.sub(this.H0.mult(deltaG));
		/*Matrix CorrectiveMatrix*/
		
		
		double CurvatureA	= deltaX.scalar(ErrorKA);
		double CurvatureH	= deltaG.scalar(ErrorKH);
		if(Math.abs(CurvatureA) >SMALL_CURVATURE && Math.abs(CurvatureH) >SMALL_CURVATURE) {
			Matrix NEWA = ErrorKA.MultEx(ErrorKA);
			this.A0	= this.A0.add(NEWA.leftmul((1/CurvatureA)));
			Matrix NEWH = ErrorKH.MultEx(ErrorKH);
			this.H0 = this.H0.add(NEWH.leftmul((1/CurvatureH)));
			System.out.println("Mise a jour des matrices effectue\n");
		}
		
		/*if(q.grad(this.iter_vec).norm() < GRADIENT_MIN_NORM || this.radius == DELTA_MIN) {
			System.out.println("STOP CONDITION FULLFILLED :");
			System.out.println("Grad val : "+f.grad(this.iter_vec).norm());
			System.out.println("Radius "+this.radius + "\n");
			throw new EndOfIteration();
		}
		**/
	}

	/* *
	* Initial trust region size .
	*/
	public final static double DELTA_INIT = 1.0;
	/* *
	* Ratio by which the region size is either
	* multiplied / divided
	*/
	public final static double DELTA_RATIO = 2.0;
	/* *
	* Minimal trust region size .
	*/
	public final static double DELTA_MIN = 1e-10;
	/* *
	* Maximal trust region size .
	*/
	public final static double DELTA_MAX = 10;
	/* *
	* Ratio of actual / predicted reduction above which
	* adequacy of the quadratic model is considered as good .
	*/
	public final static double GOOD_ADEQUACY = 0.75;
	/* *
	* Ratio of actual / predicted reduction under which
	* adequacy of the quadratic model is considered as poor .
	*/
	public final static double POOR_ADEQUACY = 0.25;
	/* *
	* Small curvature : means that the new " measurement " ( dx , dg )
	* is close to be linearly dependent of the n previous ones .
	* When the dot product between dx / dg and e is less than this value ,
	* the corresponding matrix Hk / Ak should not be updated .
	*/
	public final static double SMALL_CURVATURE = 1e-20;
	/* *
	* When the gradient norm is less than this value , the first - order
	* condition are considered to be fulfilled .
	*/
	public final static double GRADIENT_MIN_NORM = 1e-15;
}
