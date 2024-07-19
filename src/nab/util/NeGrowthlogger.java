package nab.util;

import java.io.PrintStream;
import java.util.Arrays;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;

public class NeGrowthlogger extends CalculationNode implements Loggable {
	
	
    final public Input<RealParameter> initialLogNeInput = new Input<>("initialLogNe",
            "Nes over time in log space", Input.Validate.REQUIRED);
    final public Input<RealParameter> growthInput = new Input<>("growth",
            "Nes over time in log space", Input.Validate.REQUIRED);
    final public Input<RealParameter> rateShiftsInput = new Input<>("rateShifts",
            "When to switch between elements of Ne", Input.Validate.REQUIRED);
	
    //
    // Public stuff
    //
    RealParameter initialLogNe;
    RealParameter growth;
    RealParameter rateShifts;
    
    double[] Ne;
    boolean NesKnown = false;

    @Override
	public void initAndValidate() {
    	initialLogNe = initialLogNeInput.get();    	    	
    	rateShifts = rateShiftsInput.get();
    	growth = growthInput.get();
    	Ne = new double[growth.getDimension()+1];
    }

	
	@Override
	public void init(PrintStream out) {		

		for (int i = 0 ; i < Ne.length; i++){
			out.print(String.format("Ne%d\t", i+1));
		}				
	}

	@Override
	public void log(long sample, PrintStream out) {
		recalculateNe();
		for (int i = 0 ; i < Ne.length; i++){
			out.print(Ne[i] + "\t");
		}
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}
	// computes the Ne's at the break points
	private void recalculateNe() {
		Ne = new double[growth.getDimension()+1];
		Ne[0] = initialLogNe.getValue();
		double curr_time = 0.0;
		for (int i = 1; i < Ne.length; i++) {
			Ne[i] = Ne[i-1]-growth.getArrayValue(i-1)*(rateShifts.getArrayValue(i-1)-curr_time);
			
			curr_time = rateShifts.getArrayValue(i-1);
		}
		NesKnown = true;
	}

}
