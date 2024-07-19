package nab.multitree;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import beast.core.CalculationNode;
import beast.core.Citation;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Loggable;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.coalescent.IntervalType;
import beast.mascot.distribution.MascotNative2;
import beast.mascot.dynamics.Dynamics;
import beast.mascot.ode.*;
import cern.colt.Arrays;
import nab.skygrid.TimeVaryingRates;

/**
 * @author Nicola Felix Mueller
 */

@Description("Calculates the probability of a beast.tree using under the framework of Mueller (2017).")
@Citation("Nicola F. Müller, David A. Rasmussen, Tanja Stadler (2017)\n  The Structured Coalescent and its Approximations.\n  Mol Biol Evol 2017 msx186. doi: 10.1093/molbev/msx186")
public class MultiTreeMascotLogger extends StructuredMultiTreeDistribution implements Loggable {
	
	public static boolean debug = false;
	public Input<Dynamics> dynamicsInput = new Input<>("dynamics", "Input of rates", Input.Validate.REQUIRED);
	public Input<Double> epsilonInput = new Input<>("epsilon", "step size for the RK4 integration",0.001);
	public Input<Double> maxStepInput = new Input<>("maxStep", "step size for the RK4 integration", Double.POSITIVE_INFINITY);
	
	public Input<Boolean> cacheInput = new Input<>("useCache", "use cache to speed things up (may be fragile)", false);

	enum MascotImplementation {java, indicators, allnative};
	public Input<MascotImplementation> implementationInput = new Input<>("implementation", "implementation, one of " + MascotImplementation.values().toString(),
			MascotImplementation.allnative, MascotImplementation.values());
	
    final public Input<TimeVaryingRates> immigrationRateInput = new Input<>("immigrationRate", 
    		"A population size model", Validate.REQUIRED);
    
    final public Input<Integer> minClusterSizeInput = new Input<>("minClusterSize", 
    		"A population size model", 0);


    
	public int samples;
	public int nrSamples;
    
	protected int nrLineages;   

    // current rates         
    //private double[] migrationRates;
    //private int [] indicators;
    protected double[] coalescentRates; 	

    
    // Set up for lineage state probabilities
    protected ArrayList<Integer> activeLineages;
    protected double[] linProbs;
    protected double[] linProbsNew;
    protected int linProbsLength;
    protected int states;
	
    // store the linProbs, multiplicators and logP's at coalescent points in jagged arrays from last time
    protected double[] coalLinProbs;
    protected int [] coalLinProbsLengths;
    protected double[] coalLogP;
    protected int[] coalRatesInterval;
    //private ArrayList<ArrayList<Integer>> coalActiveLineages;
    
    // deep store the things above for MCMC
    protected double[] storeLinProbs;
    protected int [] storedCoalLinProbsLengths;
    protected double[] storeLogP;
    protected int[] storeRatesInterval;
    //private ArrayList<ArrayList<Integer>> storeActiveLineages;

    
    protected double [] nextTreeEvents;
    protected double [] storedNextTreeEvents;
    protected double [] nextRateShifts;
    protected double [] storedNextRateShifts;
	//private int [] treeIntervalNrs;
	//private int [] storedTreeIntervalNrs;
	//private int [] lineagesAddded;
	//private int [] storedLineagesAddded;

    // check if this is the first calculation
    protected int first = 0;

	
	// maximum integration error tolerance
    protected double maxTolerance = 1e-3;            
    protected boolean recalculateLogP;
	Euler2ndOrderBase euler;
	public Dynamics dynamics;
	StructuredMultiTreeIntervals treeIntervals;
	
    TimeVaryingRates immigrationRate;

	Map<Integer, Integer> nodeType;

	MascotNative2 mascotImpl = null;
	boolean useCache;
	
	
	private Map<Integer, Double[]> nodeProbs;
	private Map<Integer, Double[]> rootProbs;
	
	
    @Override
    public void initAndValidate(){  	
    	useCache = cacheInput.get();
    	dynamics = dynamicsInput.get();
    	treeIntervals = multiTreeIntervalsInput.get();
    	treeIntervals.calculateIntervals();       
        nrSamples = treeIntervals.getSampleCount() + 1;    
        states = dynamics.getDimension();
                
    	int intCount = treeIntervals.getIntervalCount();

    	// initialize storing arrays and ArrayLists
    	coalLinProbs = new double[intCount * intCount * states];
    	storeLinProbs = new double[intCount * intCount * states];
    	coalLinProbsLengths = new int[intCount];
    	storedCoalLinProbsLengths = new int[intCount];
    	coalLogP = new double[intCount];
    	storeLogP = new double[intCount];
    	coalRatesInterval = new int[intCount];
    	storeRatesInterval = new int[intCount];
    	//coalActiveLineages = new ArrayList<>();
    	nextTreeEvents = new double[intCount];
    	nextRateShifts = new double[intCount];
    	storedNextTreeEvents = new double[intCount];
    	storedNextRateShifts = new double[intCount];
    	parents = new int[intCount];
    	//treeIntervalNrs = new int[intCount];
    	//storedTreeIntervalNrs = new int[intCount];
    	//lineagesAddded = new int[intCount];
    	//storedLineagesAddded = new int[intCount];
    	
    	//ArrayList<Integer> emptyList = new ArrayList<>();
    	//for (int i = 0; i <= intCount; i++) coalActiveLineages.add(emptyList);
    	
        immigrationRate = immigrationRateInput.get();

    	
    	activeLineages = new ArrayList<>();

    	int MAX_SIZE = intCount * states * 2;
    	linProbs_for_ode = new double[MAX_SIZE];
    	linProbs_tmp = new double[MAX_SIZE];
    	linProbs = new double[MAX_SIZE];
    	linProbsNew = new double[MAX_SIZE];
    	
		nodeType = new HashMap<>();
    	if (dynamics.typeTraitInput.get() != null) {
    		for (Integer number : treeIntervals.tipMapping.keySet()) {
    			nodeType.put(number, dynamics.getValue(treeIntervals.tipMapping.get(number)));
    		}
    	} else {
    		// TODO: fill in nodeType another way
    	}

    	MascotImplementation imp = implementationInput.get();
    	switch (imp) {
//    	case allnative: if (Euler2ndOrderNative.loadLibrary()) {
//    		mascotImpl = new MascotNative2(treeIntervals, nodeType, states,epsilonInput.get(), maxStepInput.get(), useCache);
//    		break;
//    	}
    	case indicators: if (Euler2ndOrderNative.loadLibrary()) {
    		euler = new Euler2ndOrderNative();
        	euler.setup(MAX_SIZE, states, epsilonInput.get(), maxStepInput.get());
        	Log.warning("Using " + euler.getClass().getSimpleName());
    		break;
    	}
    	case java:
    		switch (states) {
    		case 2: euler = new Euler2ndOrder2(); break;
    		case 3: euler = new Euler2ndOrder3(); break;
    		case 4: euler = new Euler2ndOrder4(); break;
    		case 5: euler = new Euler2ndOrder5(); break;
    		case 6: euler = new Euler2ndOrder6(); break;
    		case 7: euler = new Euler2ndOrder7(); break;
    		case 8: euler = new Euler2ndOrder8(); break;
    		case 9: euler = new Euler2ndOrder9(); break;
    		case 10: euler = new Euler2ndOrder10(); break;
    		case 11: euler = new Euler2ndOrder11(); break;
    		case 12: euler = new Euler2ndOrder12(); break;
    		case 13: euler = new Euler2ndOrder13(); break;
    		case 14: euler = new Euler2ndOrder14(); break;
    		case 15: euler = new Euler2ndOrder15(); break;
    		case 16: euler = new Euler2ndOrder16(); break;
    		case 17: euler = new Euler2ndOrder17(); break;
    		case 18: euler = new Euler2ndOrder18(); break;
    		case 19: euler = new Euler2ndOrder19(); break;
    		case 20: euler = new Euler2ndOrder20(); break;
    		case 21: euler = new Euler2ndOrder21(); break;
    		case 22: euler = new Euler2ndOrder22(); break;
    		case 23: euler = new Euler2ndOrder23(); break;
    		case 24: euler = new Euler2ndOrder24(); break;
    		case 25: euler = new Euler2ndOrder25(); break;
    		case 26: euler = new Euler2ndOrder26(); break;
    		case 27: euler = new Euler2ndOrder27(); break;
    		case 28: euler = new Euler2ndOrder28(); break;
    		case 29: euler = new Euler2ndOrder29(); break;
    		case 30: euler = new Euler2ndOrder30(); break;
    		default: euler = new Euler2ndOrder(); break;
    		}

        	euler.setup(MAX_SIZE, states, epsilonInput.get(), maxStepInput.get());
        	Log.warning("Using " + euler.getClass().getSimpleName());
    	}
    	
    	
    }
    
    double [] linProbs_for_ode;
    double [] linProbs_tmp;
    int [] parents;

    public double calculateLogP() {  	
    	
    	nodeProbs = new HashMap<>();
    	rootProbs = new HashMap<>();
    	
    	// newly calculate tree intervals (already done by swap() below)
    	treeIntervals.calculateIntervals();

        // Set up ArrayLists for the indices of active lineages and the lineage state probabilities
        activeLineages.clear();
        logP = 0;
        nrLineages = 0;
        //linProbs = new double[0];// initialize the tree and rates interval counter
        linProbsLength = 0;
        int treeInterval = 0, ratesInterval = 0;        
     	double nextEventTime = 0.0;
        
        // Time to the next rate shift or event on the tree
        double nextTreeEvent = treeIntervals.getInterval(treeInterval);
        double nextRateShift = dynamics.getInterval(ratesInterval);
        
//        if (first == 0 || !dynamics.areDynamicsKnown()) {
        	setUpDynamics();
//        }

		coalescentRates = dynamics.getCoalescentRate(ratesInterval);  
        //migrationRates = dynamics.getBackwardsMigration(ratesInterval);
		//indicators = dynamics.getIndicators(ratesInterval);
		nrLineages = activeLineages.size();
		linProbsLength = nrLineages * states;
		double currTime = 0.0;
		
		System.out.println(Arrays.toString(dynamics.getBackwardsMigration(ratesInterval)));

        // Calculate the likelihood
        do {       
        	nextEventTime = Math.min(nextTreeEvent, nextRateShift);       	
        	if (nextEventTime > 0) {													// if true, calculate the interval contribution        		
                if (recalculateLogP) {
    				System.err.println("ode calculation stuck, reducing tolerance, new tolerance= " + maxTolerance);
    				maxTolerance *=0.9;
    		    	recalculateLogP = false;
    				System.exit(0);
                	return calculateLogP();
                }
                if (nrLineages>0) {
	        		logP += doEuler(nextEventTime, ratesInterval);
	        		// add immigration term
	            	double meanMig = Math.exp(immigrationRate.getMeanRate(currTime, currTime+nextEventTime));
	            	currTime+=nextEventTime;
	            	logP -= meanMig * activeLineages.size() * nextEventTime;
                }


        	}
       	
        	if (nextTreeEvent <= nextRateShift){
 	        	if (treeIntervals.getIntervalType(treeInterval) == IntervalType.COALESCENT) {
// 	        		System.out.print(logP);

 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
	        		logP += coalesce(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);	  				// calculate the likelihood of the coalescent event
	        	}
 	       		
 	       		if (treeIntervals.getIntervalType(treeInterval) == IntervalType.SAMPLE) { 	
// 	       			System.out.print(logP);
 	       			//if (linProbsLength > 0)
 	       			//	logP += normalizeLineages(linProbs);								// normalize all lineages before event
	       			nrLineages++;													// sampling event increases the number of lineages by one
 	       			sample(treeInterval, ratesInterval, nextTreeEvent, nextRateShift);							// calculate the likelihood of the sampling event if sampling rate is given
	       		}	
 	       		
 	       		if (treeIntervals.getIntervalType(treeInterval) == IntervalType.MIGRATION) { 	
// 	       			System.out.print(logP);
 	        		nrLineages--;													// coalescent event reduces the number of lineages by one
 	        		introduction(treeInterval);	 				
                	double mig = Math.exp(immigrationRate.getRate(currTime));
                	logP += Math.log(mig);
	       		}	
 	       		
// 	       		System.out.println(logP);

 	       		
 	       		treeInterval++;
        		nextRateShift -= nextTreeEvent;   
        		try{
        			nextTreeEvent = treeIntervals.getInterval(treeInterval);
        		}catch(Exception e){
        			break;
        		}
        	} else {
        		ratesInterval++;
        		coalescentRates = dynamics.getCoalescentRate(ratesInterval);  
                //migrationRates = dynamics.getBackwardsMigration(ratesInterval);
        		//indicators = dynamics.getIndicators(ratesInterval);  
        		nextTreeEvent -= nextRateShift;
 	       		nextRateShift = dynamics.getInterval(ratesInterval);
        	}
        	if (logP == Double.NEGATIVE_INFINITY) {
        		return logP;
        	}
            if (debug) {
            	Log.info(treeInterval + " " + ratesInterval + " " + logP);
            }        	
        } while(nextTreeEvent <= Double.POSITIVE_INFINITY);
//		System.out.println();

        first++;
//        System.err.println(logP);
		return logP;  	
    }


	protected void setUpDynamics() {
    	int n = dynamics.getEpochCount();
    	double [][] coalescentRates = new double[n][];
    	double [][] migrationRates = new double[n][];
    	int [][] indicators = new int[n][];
    	double [] nextRateShift = dynamics.getIntervals();
    	for (int i = 0; i < n; i++) {
    		coalescentRates[i] = dynamics.getCoalescentRate(i);  
            migrationRates[i] = dynamics.getBackwardsMigration(i);
    		indicators[i] = dynamics.getIndicators(i);
    	}
    	dynamics.setDynamicsKnown();
		euler.setUpDynamics(coalescentRates, migrationRates, indicators, nextRateShift);

	}

	double [] storedMigrationRates = new double[0];
    double [] storedCoalescentRates = new double[0];
    int storedNrLineages = -1;
    
	public double doEuler(double nextEventTime, int ratesInterval) {
		//for (int i = 0; i < linProbs.length; i++) linProbs_tmp[i] = linProbs[i];
		if (linProbs_tmp.length != linProbsLength + 1) {
			linProbs_tmp= new double[linProbsLength + 1];
		}
		System.arraycopy(linProbs,0,linProbs_tmp,0,linProbsLength);
		linProbs_tmp[linProbsLength] = 0;

		linProbs[linProbsLength-1] = 0;
		

//		if (dynamics.hasIndicators) {
//			euler.initWithIndicators(migrationRates, indicators, coalescentRates, nrLineages);
//			euler.calculateValues(nextEventTime, linProbs_tmp, linProbsLength + 1);
//		} else {
			euler.initAndcalculateValues(ratesInterval, nrLineages, nextEventTime, linProbs_tmp, linProbsLength + 1);
//		}
		
		//		System.out.println(Arrays.toString(linProbs));		

		//for (int i = 0; i < linProbs.length; i++) linProbs[i] = linProbs_tmp[i];
		System.arraycopy(linProbs_tmp,0,linProbs,0,linProbsLength);
//		System.out.print(euler.);
		return linProbs_tmp[linProbsLength];
	}

    protected void sample(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
    	if (debug) {
    		System.err.println("sample activeLineages " + currTreeInterval + " = " + activeLineages);
    	}
		int incomingLines = treeIntervals.getLineagesAdded(currTreeInterval);
		int newLength = linProbsLength + 1 * states;
		
		int currPosition = linProbsLength;
		
		/*
		 * If there is no trait given as Input, the model will simply assume that
		 * the last value of the taxon name, the last value after a _, is an integer
		 * that gives the type of that taxon
		 */
		if (dynamics.typeTraitInput.get()!=null){
			Integer l = incomingLines; {
				activeLineages.add(l);//.getNr());
				int sampleState = nodeType.get(l);//dynamics.getValue(tree.getNode(l).getID());
				if (sampleState>= dynamics.getDimension()){
					System.err.println("sample discovered with higher state than dimension");
				}
				
				for (int i = 0; i < states; i++){
					if (i == sampleState){
						linProbs[currPosition] = 1.0;currPosition++;
					}
					else{
						linProbs[currPosition] = 0.0;currPosition++;
					}
				}
			}				
		}else{
			Integer l = incomingLines; {
//				activeLineages.add(l);//.getNr());
//				String sampleID = tree.getNode(l).getID();
				int sampleState = 0;
				if (states > 1){				
//					String[] splits = sampleID.split("_");
//					sampleState = Integer.parseInt(splits[splits.length-1]); //samples states (or priors) should eventually be specified in the XML
				}
				for (int i = 0; i < states; i++){
					if (i == sampleState){
						linProbs[currPosition] = 1.0;currPosition++;
					}
					else{
						linProbs[currPosition] = 0.0;currPosition++;
					}
				}
			}	
		}
		linProbsLength = newLength;
		// store the node
    }
          
    protected double coalesce(int currTreeInterval, int currRatesInterval, double nextTreeEvent, double nextRateShift) {
    	
    	int coalLines0 = treeIntervals.getLineagesRemoved(currTreeInterval,0);
    	int coalLines1 = treeIntervals.getLineagesRemoved(currTreeInterval,1);
    	
    	if (debug) {
    		System.err.println("coalesce activeLineages " + currTreeInterval + " " + coalLines0 + " " + coalLines1 + " = " + activeLineages);
    	}
		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines0);//.getNr());
		final int daughterIndex2 = activeLineages.indexOf(coalLines1);//.getNr());
		
		if (daughterIndex1 == -1 || daughterIndex2 == -1) {
			System.out.println(coalLines0/*.getNr()*/ + " " + coalLines1/*.getNr()*/ + " " + activeLineages);
			System.out.println("daughter lineages at coalescent event not found");
			System.exit(0);
			return Double.NaN;
		}
		
		
		double[] lambda = new double[states];
		double lambdaSum = 0.0;
		
		/*
		 * Calculate the overall probability for two strains to coalesce 
		 * independent of the state at which this coalescent event is 
		 * supposed to happen
		 */
        for (int k = 0; k < states; k++) { 
        	Double pairCoalRate = coalescentRates[k] * linProbs[daughterIndex1*states + k] * linProbs[daughterIndex2*states + k];			
			if (!Double.isNaN(pairCoalRate)){
				lambda[k] = pairCoalRate;
				lambdaSum+=pairCoalRate;
			} else {
				return Double.NEGATIVE_INFINITY;
			}
        }

        int lineageToAdd = treeIntervals.getLineagesAdded(currTreeInterval);
        activeLineages.add(lineageToAdd);        
		
		int linCount = 0;		
		// add all lineages execpt the daughter lineage to the new p array
		for (int i = 0; i <= nrLineages; i++){
			if (i != daughterIndex1 && i != daughterIndex2){
				for (int j = 0; j < states; j++){
					linProbsNew[linCount*states + j] = linProbs[i*states + j];
				}
				linCount++;
			}
		}
		// add the parent lineage
		Double[] nodeProbs = new Double[lambda.length];
		for (int j = 0; j < states; j++){
			linProbsNew[linCount*states + j] = lambda[j]/lambdaSum;
			nodeProbs[j] = lambda[j]/lambdaSum;
		}
		
		this.nodeProbs.put(lineageToAdd, nodeProbs);
		
		
		// set p to pnew
		//double [] tmp = linProbs;
		linProbs = linProbsNew;	
		linProbsNew = linProbs;
		linProbsLength = linProbsLength - states;
		
		
		//Remove daughter lineages from the line state probs
		if (daughterIndex1>daughterIndex2){
			// remove the daughter lineages from the active lineages
			activeLineages.remove(daughterIndex1);
			activeLineages.remove(daughterIndex2);			
		} else {
			// remove the daughter lineages from the active lineages
			activeLineages.remove(daughterIndex2);
			activeLineages.remove(daughterIndex1);			
		}		
     
		
		if (lambdaSum==0)
			return Double.NEGATIVE_INFINITY;
		else
			return Math.log(lambdaSum);
    }
    
	private void introduction(int treeInterval) {
    	int coalLines0 = treeIntervals.getLineagesRemoved(treeInterval, 0);	 	   		
    	final int daughterIndex1 = activeLineages.indexOf(coalLines0);//.getNr());
    	int linCount=0;
    	Double[] rootProbs = new Double[states];
   		for (int i = 0; i <= nrLineages; i++){
			if (i != daughterIndex1){
				for (int j = 0; j < states; j++){
					linProbsNew[linCount*states + j] = linProbs[i*states + j];
				}
				linCount++;
			}else {
				for (int j = 0; j < states; j++){
					rootProbs[j] = linProbs[i*states + j];
				}
			}
		}
   		
   		
   		this.rootProbs.put(daughterIndex1, rootProbs);
   		
		// set p to pnew
		//double [] tmp = linProbs;
		linProbs = linProbsNew;	
		linProbsNew = linProbs;
		linProbsLength = linProbsLength - states;
    	
    	
		activeLineages.remove(daughterIndex1);		
	}

       
    /**
     * Loggable interface implementation follows *
     */
    @Override
    public void init(final PrintStream out) {
        out.println("#NEXUS\n");
        out.println("Begin trees;\n");        
    }

    @Override
    public void log(final long sample, final PrintStream out) {
        calculateLogP();
        double maxHeight=-1.0;
        for (int i = 0; i < treeIntervals.treeInput.get().size();i++) {
        	maxHeight = Math.max(treeIntervals.treeInput.get().get(i).getRoot().getHeight()+
        			treeIntervals.offset[i] + treeIntervals.rootLengthInput.get().get(i).getValue(), 
        			maxHeight);
        }
        String tree_string = "rem";
        
        for (int i = 0; i < treeIntervals.treeInput.get().size();i++) {
        	if (treeIntervals.treeInput.get().get(i).getExternalNodes().size()>=minClusterSizeInput.get())
	        	tree_string = tree_string + "," + 
		        getNewick(treeIntervals.treeInput.get().get(i), treeIntervals.treeNodeCount[i], i) + 
		        ":" +  (maxHeight-treeIntervals.treeInput.get().get(i).getRoot().getHeight()-treeIntervals.offset[i]-treeIntervals.rootLengthInput.get().get(i).getValue());
        }
        
        tree_string = tree_string.replace("rem,", "(");
        tree_string = tree_string + "):0.0";
        
        
//        System.out.println(Arrays.toString(nodeProbs.get(96)));
       
        
        out.print("tree STATE_" + sample + " = ");
        out.print(tree_string);
        out.print(";");	        		
    }

    private String getNewick(Tree tree, int nodeOffset, int i) {
        final StringBuilder rootString = new StringBuilder();

        rootString.append("(" + toNewick(tree.getRoot(), nodeOffset));
        rootString.append(":" + treeIntervals.rootLengthInput.get().get(i).getValue());
        rootString.append(")[&obs=0]");
//        rootString.append(")" + getStateProbs(rootProbs.get(tree.getRoot().getNr()+nodeOffset)));
		return rootString.toString();
	}
    
    public String toNewick(Node n, int nodeOffset) {
        final StringBuilder buf = new StringBuilder();
        if (!n.isLeaf()) {
            buf.append("(");
            boolean isFirst = true;
            for (Node child : n.getChildren()) {
                if (isFirst)
                    isFirst = false;
                else
                    buf.append(",");
                
                buf.append(toNewick(child, nodeOffset));
            }
            buf.append(")");
            buf.append(getStateProbs(nodeProbs.get(n.getNr()+nodeOffset)));            

            if (n.getID() != null)
                buf.append(n.getID());
        } else {
            buf.append(n.getID());
            buf.append(getLeafProbs(nodeType.get(n.getNr()+nodeOffset)));            

        }
        if (!n.isRoot())
        	buf.append(":").append(n.getLength());
        return buf.toString();
    }
    
    private String getLeafProbs(int nodeType) {
    	final StringBuilder meta =  new StringBuilder();
    	meta.append("[&obs=1,max=");
    	
        meta.append(dynamicsInput.get().getStringStateValue(nodeType));
        for (int i = 0; i < states; i++) {
        	meta.append("," + dynamicsInput.get().getStringStateValue(i));
        	if (i==nodeType) {
	        	meta.append("=1.0");
        	}else {
        		meta.append("=0.0");
        	}
        }
        meta.append("]");
        return meta.toString();
	}


	private String getStateProbs(Double[] probs) {
    	final StringBuilder meta =  new StringBuilder();
    	meta.append("[&obs=1,max=");
    	
        double maxVal =-1;
        int maxState =-1;
        for (int i = 0; i < states; i++) {
        	if (probs[i]>maxVal) {
        		maxVal = probs[i];
        		maxState=i;        				
        	}
        }
        meta.append(dynamicsInput.get().getStringStateValue(maxState));
        for (int i = 0; i < states; i++) {
        	meta.append("," + dynamicsInput.get().getStringStateValue(i));
        	meta.append("=" + probs[i]);
        }
        meta.append("]");
        return meta.toString();
    } 


	@Override
    public void close(final PrintStream out) {
        out.print("End;");
    }

}
