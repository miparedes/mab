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

public class ClusterClock extends BranchRateModel.Base implements Loggable {
	
    public final Input<Tree> treeInput = new Input<>("tree", "the tree containing the taxon set", Validate.REQUIRED);
	final public Input<RealParameter> meanClockRateInput = new Input<>("meanClockRate",
			"the mean Clock Rate.", Input.Validate.REQUIRED);
	final public Input<RealParameter> relativeClockRateInput = new Input<>("relativeClockRate",
			"the relative Clock Rate of an individual cluster.", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> rateCategoryInput = new Input<>("rateCategory",
			"the rate category a cluster belongs to.", Input.Validate.REQUIRED);	
    public final Input<List<TaxonSet>> taxonsetInput = new Input<>("taxonset",
            "set of taxa set to define an individual clusters", new ArrayList());
	final public Input<BooleanParameter> useRootBranchInput = new Input<>("useRootBranch",
			"if true, the rate is for the branch above the root of the taxon set mrca.", Input.Validate.REQUIRED);
	final public Input<IntegerParameter> rateMapInput = new Input<>("rateMap",
			"maps taxon set to rate category to allow multiple sets to a priori have the same rate.", Input.Validate.REQUIRED);

    
    int nr_categories;
    RealParameter mean;
    RealParameter relative;
    IntegerParameter category;
    Tree tree;
    
    // array of indices of taxa
    int[][] taxonIndex;
    
	List<Integer> nodevals;
	List<Integer> storedNodevals;
	
	List<Integer> cluster;
	List<Integer> storedCluster;
	
	int[] rateMap;
	
    double[] unscaledBranchRates;
    double scaleFactor = 1.0;


	
	@Override
	public void initAndValidate() {
		rateMap = new int[rateMapInput.get().getDimension()];
		int max_val = -1;
		for (int i = 0; i < rateMap.length; i++) {
			rateMap[i] = (int) rateMapInput.get().getArrayValue(i);
			max_val = Math.max(max_val, rateMap[i]);
		}
			
		
		nr_categories = max_val+2;
		mean = meanClockRateInput.get();
		relative = relativeClockRateInput.get();
		category = rateCategoryInput.get();
		tree = treeInput.get();
		
		relative.setDimension(nr_categories);
		category.setDimension(nr_categories);
		category.setLower(0);
		category.setUpper(nr_categories-1);
		
		taxonIndex = new int[rateMap.length][];		
		
        final List<String> taxaNames = new ArrayList<>();
        
        for (final String taxon : tree.getTaxaNames()) {
            taxaNames.add(taxon);
        }
		
		for (int i = 0; i < taxonsetInput.get().size(); i++) {
	        List<String> set = null;
            set = taxonsetInput.get().get(i).asStringList();
            taxonIndex[i] = new int[set.size()];
            int k = 0;
	        for (final String taxon : set) {
	            final int taxonIndex_ = taxaNames.indexOf(taxon);
	            if (taxonIndex_ < 0) {
	                throw new RuntimeException("Cannot find taxon " + taxon + " in data");
	            }
	            taxonIndex[i][k++] = taxonIndex_;
	        }
		}
		
        unscaledBranchRates = new double[tree.getNodeCount()];

		recomputeClusterMemberships();	
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
		recomputeClusterMemberships();
		recalculateScaleFactor();
		return true;
	}

	private void recomputeClusterMemberships() {
		nodevals = new ArrayList<>();
		cluster = new ArrayList<>();
		for (int i = 0; i < taxonIndex.length; i++) {
			if (useRootBranchInput.get().getArrayValue(i)>0.5) {
				Node mrca = getCommonAncestorInternal(i, -1);
				nodevals.add(mrca.getNr());
				cluster.add(rateMap[i]);
			}else {
				getCommonAncestorInternal(i, rateMap[i]);
			}					
		}
	}
	
    private Node getCommonAncestorInternal(int index, int index2) {
    	Node cur = tree.getNode(taxonIndex[index][0]);

        for (int k = 1; k < taxonIndex[index].length; ++k) {
            cur = getCommonAncestor(cur, tree.getNode(taxonIndex[index][k]), index2);
        }
    
        return cur;
    }
    
    private void calculateUnscaledBranchRates(Node node) {
    	
		int index = nodevals.indexOf(node.getNr());
		double rate;
		if (index==-1) {// return background rate
			rate = relative.getArrayValue(category.getValue(nr_categories-1));			
		}else { // return the rate of the rate categorys
			rate = relative.getArrayValue(category.getValue(cluster.get(index)));			
		}


        unscaledBranchRates[node.getNr()] = rate;

        if (!node.isLeaf()) {
            calculateUnscaledBranchRates(node.getLeft());
            calculateUnscaledBranchRates(node.getRight());
        }
    }

    
    private void recalculateScaleFactor() {
        calculateUnscaledBranchRates(tree.getRoot());


        double timeTotal = 0.0;
        double branchTotal = 0.0;

        for (int i = 0; i < tree.getNodeCount(); i++) {
            Node node = tree.getNode(i);
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


    protected Node getCommonAncestor(Node n1, Node n2, int index) {
        while (n1 != n2) {
        	if (index!=-1) {
	        	if (!nodevals.contains(n1.getNr())) {
	        		nodevals.add(n1.getNr());
	        		cluster.add(index);
	        	}
	        	if (!nodevals.contains(n2.getNr())) {
	        		nodevals.add(n2.getNr());
	        		cluster.add(index);
	        	}
        	}

        		
	        double h1 = n1.getHeight();
	        double h2 = n2.getHeight();
	        if ( h1 < h2 ) {
	            n1 = n1.getParent();
	        } else if( h2 < h1 ) {
	            n2 = n2.getParent();
	        } else {
	            //zero length branches hell
	            Node n;
	            double b1 = n1.getLength();
	            double b2 = n2.getLength();
	            if( b1 > 0 ) {
	                n = n2;
	            } else { // b1 == 0
	                if( b2 > 0 ) {
	                    n = n1;
	                } else {
	                    // both 0
	                    n = n1;
	                    while( n != null && n != n2 ) {
	                        n = n.getParent();
	                    }
	                    if( n == n2 ) {
	                        // n2 is an ancestor of n1
	                        n = n1;
	                    } else {
	                        // always safe to advance n2
	                        n = n2;
	                    }
	                }
	            }
	            if( n == n1 ) {
                    n = n1 = n.getParent();
                } else {
                    n = n2 = n.getParent();
                }
	        }
        }
        return n1;
    }
    
    @Override
    public void restore() {
    	recomputeClusterMemberships();
    	recalculateScaleFactor();
        super.restore();
    }


	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < taxonsetInput.get().size(); i++) {
			out.print(String.format("rate.%s\t", taxonsetInput.get().get(i).getID()));
		}
		out.print(String.format("rate.other\t"));

	}

	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0; i < taxonsetInput.get().size(); i++) {
			double rate = relative.getArrayValue(category.getValue(rateMap[i])) * scaleFactor;
			out.print(rate + "\t");
		}
		double rate = relative.getArrayValue(category.getValue(nr_categories-1)) * scaleFactor;
		out.print(rate + "\t");
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
	}

}
