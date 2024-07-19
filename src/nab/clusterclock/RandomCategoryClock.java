package nab.clusterclock;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import cern.colt.Arrays;

public class RandomCategoryClock extends BranchRateModel.Base {
	
    public final Input<Tree> treeInput = new Input<>("tree", "the tree containing the taxon set", Validate.REQUIRED);
	final public Input<RealParameter> meanClockRateInput = new Input<>("meanClockRate",
			"the mean Clock Rate.", Input.Validate.REQUIRED);
	final public Input<RealParameter> relativeClockRateInput = new Input<>("relativeClockRate",
			"the relative Clock Rate of an individual cluster.", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> rateCategoryInput = new Input<>("rateCategory",
			"the rate category a cluster belongs to.", Input.Validate.REQUIRED);	

    
    int nr_categories;
    RealParameter mean;
    RealParameter relative;
    IntegerParameter category;
       	
    double[] unscaledBranchRates;
    double scaleFactor = 1.0;


	
	@Override
	public void initAndValidate() {		
		
		mean = meanClockRateInput.get();
		relative = relativeClockRateInput.get();
		category = rateCategoryInput.get();
		
		category.setDimension(treeInput.get().getNodeCount());
		category.setLower(0);
		category.setUpper(relative.getDimension()-1);		
		
        unscaledBranchRates = new double[treeInput.get().getNodeCount()];

		recalculateScaleFactor();
	}

	@Override
	public double getRateForBranch(Node node) {	
		if (node.isRoot()) {
			return 0.0;		
		}	
        return unscaledBranchRates[node.getNr()] * scaleFactor;
	}

	@Override
	protected boolean requiresRecalculation() {
		recalculateScaleFactor();
		return true;
	}

    
    private void calculateUnscaledBranchRates(Node node) {
        unscaledBranchRates[node.getNr()] = relative.getArrayValue(category.getValue(node.getNr()));

        if (!node.isLeaf()) {
            calculateUnscaledBranchRates(node.getLeft());
            calculateUnscaledBranchRates(node.getRight());
        }
    }

    
    private void recalculateScaleFactor() {
        calculateUnscaledBranchRates(treeInput.get().getRoot());


        double timeTotal = 0.0;
        double branchTotal = 0.0;

        for (int i = 0; i < treeInput.get().getNodeCount(); i++) {
            Node node = treeInput.get().getNode(i);
            if (!node.isRoot()) {

                double branchInTime = node.getParent().getHeight() - node.getHeight();

                double branchLength = branchInTime * unscaledBranchRates[node.getNr()];

                timeTotal += branchInTime;
                branchTotal += branchLength;
            }
        }

        scaleFactor = timeTotal / branchTotal;

        scaleFactor *= mean.getValue();
    }

    
    @Override
    public void restore() {
    	recalculateScaleFactor();
        super.restore();
    }




}
