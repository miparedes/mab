package nab.skygrid;

import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.coalescent.PopulationFunction;


/**
 * @author Nicola F. Mueller
 */
@Description("Skygrid style population prior that .")
public class MultiEpochSkygrid extends PopulationFunction.Abstract {
	
    final public Input<RealParameter> logNeInput = new Input<>("logNe",
            "Nes over time in log space", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts",
            "When to switch between elements of Ne", Input.Validate.REQUIRED);

    RealParameter logNe;
    RealParameter rateShifts;
    

    @Override
	public void initAndValidate() {
    	logNe = logNeInput.get();
    	rateShifts = rateShiftsInput.get();
    	logNe.setDimension(rateShifts.getDimension());    	
    }


	@Override
	public List<String> getParameterIds() {
		throw new IllegalArgumentException("Not implemented");
	}


	@Override
	public double getPopSize(double t) {				
//		return logNe.getArrayValue(getIntervalNr(t));
		int intervalnr = getIntervalNr(t);
		if (intervalnr>=rateShifts.getDimension()) {
//			System.err.println("index out of range, return neg inf");
			return Double.POSITIVE_INFINITY;
		}
		

		return Math.exp(logNe.getArrayValue(intervalnr));
	}

	
    @Override
	public double getIntegral(double start, double finish) {
    	if (start==finish)
    		return 0.0;
    	
    	// get the interval "start" is in
    	int first_int = getIntervalNr(start);
    	// get the interval "finish" is in
    	int last_int = getIntervalNr(finish);
    	
    	double weighted = 0.0;
    	double curr_time = start; 
    	
//    	System.out.println(first_int + " " + last_int + " " + finish + " " + start);

    	for (int i = first_int; i <= last_int;i++) {
    		if (i>=rateShifts.getDimension()) {
//    			System.err.println("index out of range, return neg inf");
    			return Double.POSITIVE_INFINITY;
    		}
    		
    		double next_time = Math.min(getTime(i), finish);
    		weighted +=  (next_time - curr_time)/ Math.exp(logNe.getArrayValue(i));
    		curr_time = next_time;
    	}
    	return weighted;
    }



	@Override
	public double getInverseIntensity(double x) {
		throw new IllegalArgumentException("Not implemented");
	}
	
	private int getIntervalNr(double t) {
		// check which interval t + offset is in
		for (int i = 0; i < rateShifts.getDimension(); i++)
			if (t<rateShifts.getArrayValue(i))
				return i;
		
		// after the last interval, just keep using the last element
		return rateShifts.getDimension();					
	}

	private double getTime(int i) {
		return rateShifts.getArrayValue(i);
	}



	@Override
	public double getIntensity(double t) {
		throw new IllegalArgumentException("Not implemented");
	}

}