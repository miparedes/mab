package nab.clusterclock;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.util.Log;
import beast.util.Randomizer;


@Description("A uniform random operator that selects a random dimension of the integer parameter and picks a new random value within the bounds.")
public class IntUniformOperator extends Operator {
    final public Input<IntegerParameter> parameterInput = new Input<>("parameter", "the parameter to operate a random walk on.", Validate.REQUIRED);
    final public Input<BooleanParameter> indicatorInput = new Input<>("indicator", "indicates which of the dimension " +
            "of the parameters can be scaled. Only used when scaleAllIndependently=false and scaleAll=false. If not specified " +
            "it is assumed all dimensions are allowed to be scaled.");


    @Override
	public void initAndValidate() {    	
        final BooleanParameter indicators = indicatorInput.get();
        if (indicators != null) {
            final int dataDim = parameterInput.get().getDimension();
            final int indsDim = indicators.getDimension();
            if (!(indsDim == dataDim || indsDim + 1 == dataDim)) {
                throw new IllegalArgumentException("indicator dimension not compatible from parameter dimension");
            }
        }

    	
    	
    }

    /**
     * override this for proposals,
     * returns log of hastingRatio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
    	
        final int index;
        final BooleanParameter indicators = indicatorInput.get();
        IntegerParameter param = parameterInput.get(this);

        final int dim = param.getDimension();

    	if (indicators != null) {
            final int dimCount = indicators.getDimension();
            final Boolean[] indicator = indicators.getValues();
            final boolean impliedOne = dimCount == (dim - 1);

            // available bit locations. there can be hundreds of them. scan list only once.
            final int[] loc = new int[dimCount + 1];
            int locIndex = 0;

            if (impliedOne) {
                loc[locIndex] = 0;
                ++locIndex;
            }
            for (int i = 0; i < dimCount; i++) {
                if (indicator[i]) {
                    loc[locIndex] = i + (impliedOne ? 1 : 0);
                    ++locIndex;
                }
            }

            if (locIndex > 0) {
                final int rand = Randomizer.nextInt(locIndex);
                index = loc[rand];
            } else {
                return Double.NEGATIVE_INFINITY; // no active indicators
            }

        } else {
            // any is good
            index = Randomizer.nextInt(dim);
        }
    	
        int newValue = Randomizer.nextInt(param.getUpper() - param.getLower() + 1) + param.getLower();
        param.setValue(index, newValue);

        return 0.0;
    }

    @Override
    public void optimize(double logAlpha) {
        // nothing to optimise
    }

} // class IntUniformOperator