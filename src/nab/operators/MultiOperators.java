package nab.operators;


import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;


@Description("Scales a parameter or a complete beast.tree (depending on which of the two is specified.")
public class MultiOperators extends Operator {
    final public Input<List<Operator>> operatorsInput = new Input<>("operator",
            "multiple operators to run at the same time", new ArrayList<>());
    
    final public Input<RealParameter> rep = new Input<>("parameter",
            "multiple operators to run at the same time",Input.Validate.REQUIRED);


    final public Input<String> scale = new Input<>("scaleFactor",
            "multiple operators to run at the same time", "");

    
    Double[][] values; 

    @Override
    public void initAndValidate() {
    	
    }



    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
    	double logHR = 0.0;
    	for (Operator o : operatorsInput.get())
    		logHR += o.proposal();   			
    	
    	return logHR;
    }
    
    @SuppressWarnings("unchecked")
	protected void readLogFile(String fileName, int burnInPercentage, String[] useLabels) throws IOException {

    } // readLogFile



} // class ScaleOperator
