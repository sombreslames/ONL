package func;

import util.Vector;

public class Rosenbrock implements RealFunc {


	public double eval(Vector v) {
		
		/* TODO */
		double x = v.get(0);
		double y = v.get(1);
		return (Math.pow((1-x),2) + 100*Math.pow((y-Math.pow(x, 2)),2));
	}

	
	public Vector grad(Vector v) {
		double x = v.get(0);
		double y = v.get(1);
		/*v[0] -- > x v[1] --> y */ 
		Vector Grad = new Vector(v.size());
		/*df/dx -2 +2x - 400x(y-x^2)*/
		Grad.set(0,-2*(1-x) -(400*x)*(y-Math.pow(x, 2)));
		/* df/dx 200*(y-x^2) */
		Grad.set(1,200*(y-Math.pow(x, 2)));
		return Grad;
	}
	
	@Override
	public int dim() {
		return 2;
	}

	
}