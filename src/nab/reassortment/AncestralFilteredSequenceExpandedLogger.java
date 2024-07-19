package nab.reassortment;



import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Input.Validate;
import beast.evolution.datatype.DataType;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.FilteredAlignment;

@Description("Reconstructs sequences at internal nodes and logs them in NEXUS format")
public class AncestralFilteredSequenceExpandedLogger extends BEASTObject implements Function, Loggable {
	
    public Input<String> tagInput = new Input<String>("tag","label used to report trait", Validate.REQUIRED);

	final public Input<Boolean> logMutationsOnlyInput = new Input<>("logMutationsOnly", "if true, logs only mutations on edges instead of sequences on nodes", false);
	
	final public Input<List<AncestralStateExpandedTreeLikelihood>> ancestralTreeLikelihoodInput = new Input<>("ancestralTreeLikelihood", "ancestrsal tree likelihoods that," +
	"if more than one is inputted is assumed to be a filtered alignment", new ArrayList<>());
	
	final public Input<Boolean> translateInput = new Input<>("translate", "if true, dna is translated to amino acids", false);

	final public Input<Boolean> relativeToRootInput = new Input<>("relativeToRoot", "if true, difference is computed relative to root", false);

	final public Input<List<String>> readingFrameInput = new Input<>("readingFrame", "denotes different reading frames", new ArrayList<>());

	final public Input<Boolean> computeTOSInput = new Input<>("computeTOS", "if true, computes the time of survival for each node into the future", false);
	
	final public Input<Boolean> logFullSequenceInput = new Input<>("logFullSequence", "if true, logs the full sequence at every node", false);
	
	final public Input<Boolean> forTreeioInput = new Input<>("forTreeio", "if true, logs the mutations to be read in by treeio", false);

	
//	int [] siteStates;
	
    protected DataType dataType;
    
    int totalLength;
    List<Integer[]> mapping;
    
    Map<String, String> translationsTable;

    int[] from;
    int[] to;
    
	@Override
	public void initAndValidate() {
		if (ancestralTreeLikelihoodInput.get().get(0).dataInput.get() instanceof FilteredAlignment) {
			totalLength = 0;
		}else {
			totalLength = ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getSiteCount();
		}
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			if (ancestralTreeLikelihoodInput.get().get(i).dataInput.get() instanceof FilteredAlignment) {
				totalLength = ((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput.get()).alignmentInput.get().getSiteCount();			
			}
		}

		if (readingFrameInput.get().size()>0) {
			from = new int[readingFrameInput.get().size()];
			to = new int[readingFrameInput.get().size()];
			for (int i = 0; i < from.length; i++) {
				String[] tmp = readingFrameInput.get().get(i).split("\\s+");
				from[i] = Integer.parseInt(tmp[0]);
				to[i] = Integer.parseInt(tmp[1]);
			}			
		}else {
			from = new int[1];
			to = new int[1];
			from[0] = 0;
			to[0] = totalLength;
		}
		
	}	

	@Override
	public void init(PrintStream out) {
		ancestralTreeLikelihoodInput.get().get(0).updateTree();
		ancestralTreeLikelihoodInput.get().get(0).initialize();
		((Tree) ancestralTreeLikelihoodInput.get().get(0).tree).init(out);
		
		dataType =  ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getDataType();

		mapping = new ArrayList<>();
		
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			if (ancestralTreeLikelihoodInput.get().get(i).dataInput.get() instanceof FilteredAlignment) {
				mapping.add(parseFilter(((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput.get()).filterInput.get(),
						((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput.get()).alignmentInput.get(),
						ancestralTreeLikelihoodInput.get().get(i).dataInput.get().getSiteCount()));
			}else {
				Integer[] newMap = new Integer[totalLength];
				for (int j = 0; j < newMap.length; j++)
					newMap[j] = j;
				mapping.add(newMap);
			}
		}
		if (translateInput.get())
			initTable();		
	}
	
	
    private Integer[] parseFilter(String filterString, Alignment alignment, int redCount) {
        // parse filter specification
        String[] filters = filterString.split(",");
        int[] from = new int[filters.length];
        int[] to = new int[filters.length];
        int[] step = new int[filters.length];
        for (int i = 0; i < filters.length; i++) {
            filterString = " " + filters[i] + " ";
            if (filterString.matches(".*-.*")) {
                // range, e.g. 1-100/3
                if (filterString.indexOf('\\') >= 0) {
                	String str2 = filterString.substring(filterString.indexOf('\\') + 1); 
                	step[i] = parseInt(str2, 1);
                	filterString = filterString.substring(0, filterString.indexOf('\\'));
                } else {
                	step[i] = 1;
                }
                String[] strs = filterString.split("-");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], alignment.getSiteCount()) - 1;
            } else if (filterString.matches(".*:.*:.+")) {
                // iterator, e.g. 1:100:3
                String[] strs = filterString.split(":");
                from[i] = parseInt(strs[0], 1) - 1;
                to[i] = parseInt(strs[1], alignment.getSiteCount()) - 1;
                step[i] = parseInt(strs[2], 1);
            } else if (filterString.trim().matches("[0-9]*")) {
                from[i] = parseInt(filterString.trim(), 1) - 1;
                to[i] = from[i];
            	step[i] = 1;
            } else {
                throw new IllegalArgumentException("Don't know how to parse filter " + filterString);
            }
        }

        boolean[] used = new boolean[alignment.getSiteCount()];
        for (int i = 0; i < to.length; i++) {
            for (int k = from[i]; k <= to[i]; k += step[i]) {
            	used[k] = true;
            }
        }
        
        Integer[] map = new Integer[redCount];
        int c=0;
        for (int i = 0; i < used.length; i++) {
        	if (used[i]) {
        		map[c] = i;
        		c++;
        	}
        }

        return map;
    }
    
    int parseInt(String str, int defaultValue) {
        str = str.replaceAll("\\s+", "");
        try {
            return Integer.parseInt(str);
        } catch (Exception e) {
            return defaultValue;
        }
    }


	boolean stop = false;
	@Override
	public void log(long sample, PrintStream out) {
		
		dataType =  ancestralTreeLikelihoodInput.get().get(0).dataInput.get().getDataType();

		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			if (i==0) {
				ancestralTreeLikelihoodInput.get().get(i).updateTree();
			}else {
				ancestralTreeLikelihoodInput.get().get(i).setTree(ancestralTreeLikelihoodInput.get().get(0).tree);
			}
			ancestralTreeLikelihoodInput.get().get(i).calculateLogP();
		}
		
		mapping = new ArrayList<>();

		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
			if (ancestralTreeLikelihoodInput.get().get(i).dataInput.get() instanceof FilteredAlignment) {
				mapping.add(parseFilter(((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput.get()).filterInput.get(),
						((FilteredAlignment) ancestralTreeLikelihoodInput.get().get(i).dataInput.get()).alignmentInput.get(),
						ancestralTreeLikelihoodInput.get().get(i).dataInput.get().getSiteCount()));
			}else {
				Integer[] newMap = new Integer[totalLength];
				for (int j = 0; j < newMap.length; j++)
					newMap[j] = j;
				mapping.add(newMap);
			}
		}
		
		if (translateInput.get())
			initTable();		

		
        out.print("tree STATE_" + sample + " = ");
        Tree tree = ancestralTreeLikelihoodInput.get().get(0).tree;
        tree.getRoot().sort();
        
        if (computeTOSInput.get())
        	getTOS(tree);

        out.print(toNewick(tree.getRoot(), null, null));
        out.print(";");
        if (stop) {
        	System.out.println(ancestralTreeLikelihoodInput.get().get(0).networkInput.get());
        	System.out.println(toNewick(tree.getRoot(), null, null));
        	System.exit(0);
        }
	}

	
    String toNewick(Node node, int[] parentSiteStates, List<List<String>> translatedParent) {
        StringBuffer buf = new StringBuffer();
        
        int[] currSiteStates = new int[totalLength];
        List<List<String>> translated = new ArrayList<>();
        
        for (int i = 0; i < currSiteStates.length; i++)
        	currSiteStates[i] =-1;
        
		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
		    int [] patternstates = ancestralTreeLikelihoodInput.get().get(i).getStatesForNode(
		    		ancestralTreeLikelihoodInput.get().get(0).tree, node);
		    
		    for (int j = 0; j < mapping.get(i).length; j++) {
		    	currSiteStates[mapping.get(i)[j]] = patternstates[j];
		    }	    
		}
		
	    if (translateInput.get()) {
	    	translate(currSiteStates, translated);
	    	if (translatedParent==null)
	    		translatedParent = new ArrayList<>(translated);
	    }
	    if (relativeToRootInput.get()) {
	    	if (parentSiteStates==null) {
	    		parentSiteStates = new int[currSiteStates.length];
	    		System.arraycopy(currSiteStates, 0, parentSiteStates, 0, parentSiteStates.length);
	    	}
	    }

        if (node.getLeft() != null) {
            buf.append("(");
            if (relativeToRootInput.get())
            	buf.append(toNewick(node.getLeft(), parentSiteStates, translatedParent));
            else 
            	buf.append(toNewick(node.getLeft(), currSiteStates, translated));
            if (node.getRight() != null) {
                buf.append(',');
                if (relativeToRootInput.get())
                	buf.append(toNewick(node.getRight(), parentSiteStates, translatedParent));
                else
                	buf.append(toNewick(node.getRight(), currSiteStates, translated));
            }
            
            buf.append(")");
        } else {
            buf.append(node.getNr() + 1);
        }
        
//		
//		for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
//		    int [] patternstates = getStatesForNode(treeInput.get(), node);
//		    for (int j = 0; j < siteStates.length; j++) {
//		    	siteStates[j] = patternstates[dataInput.get().getPatternIndex(i)];
//		    }

	    if (logMutationsOnlyInput.get()) {
	    	if (parentSiteStates!=null) {

			    // look for differences
	    		if (!relativeToRootInput.get()) {
	    			if (!translateInput.get()) {
					    List<Integer> diffs = new ArrayList<>();
					    for (int i = 0; i < parentSiteStates.length; i++) {
					    	if (parentSiteStates[i]!=currSiteStates[i]) {
					    		diffs.add(i);
					    	}
					    	
					    }
					    if (diffs.size()>0) {
					    	buf.append("[&" + tagInput.get() + "=\"");
						    for (int i = 0; i < diffs.size(); i++) {
						    	buf.append(dataType.getCharacter(parentSiteStates[diffs.get(i)]) + "" + (diffs.get(i)+1) + 
						    			"" + dataType.getCharacter(currSiteStates[diffs.get(i)]));
						    	if (i<(diffs.size()-1)) {
						    		buf.append(",");
						    	}
						    }  
						    buf.append("\"");
					    	buf.append(",mutations="+diffs.size());
					    	if (node.getMetaData("distance")!=null)
						    	buf.append(",tos="+node.getMetaData("distance"));
						    buf.append("]");
					    }
	    			}else {
					    List<Integer> diffs = new ArrayList<>();
					    for (int i = 0; i < translatedParent.size(); i++) {
					    	if (translatedParent.get(i)!=translated.get(i)) {
					    		diffs.add(i);
					    	}
					    	
					    }
					    if (diffs.size()>0) {
					    	buf.append("[&" + tagInput.get() + "=\"");
						    for (int i = 0; i < diffs.size(); i++) {
						    	buf.append(translatedParent.get(diffs.get(i)) + "" + (diffs.get(i)+1) + 
						    			"" + translated.get(diffs.get(i)));
						    	if (i<(diffs.size()-1)) {
						    		buf.append(",");
						    	}
						    }  
						    buf.append("\"");
					    	buf.append(",mutations="+diffs.size());
					    	if (node.getMetaData("distance")!=null)
						    	buf.append(",tos="+node.getMetaData("distance"));

						    buf.append("]");
					    }

	    			}
	    		}else {
				    List<Integer> diffs = new ArrayList<>();
				    for (int i = 0; i < parentSiteStates.length; i++) {
				    	if (parentSiteStates[i]!=currSiteStates[i]) {
				    		diffs.add(i);
				    	}
				    	
				    }
				    	
			    	buf.append("[&NT_muts=" + diffs.size());
	    			if (translateInput.get()) {
					    diffs = new ArrayList<>();
					    int c=0;
					    for (int i = 0; i < translatedParent.size(); i++) {
					    	for (int j = 0; j < translatedParent.get(i).size(); j++) {
						    	if (translatedParent.get(i).get(j)!=translated.get(i).get(j)) {
						    		c++;
						    	}
					    	}
					    }
				    	buf.append(",AA_muts=" + c);
	    			}
			    	if (node.getMetaData("distance")!=null)
				    	buf.append(",tos="+node.getMetaData("distance"));
			    	
			    	if (node.isLeaf()) {
				    	if (node.getID().startsWith("reaSurvived_")) {
					    	buf.append(",co={");
					    	String co_seg = "rem";
					    	for (int i = 0; i < ((BitSet) node.getMetaData("co")).size(); i++) {
					    		if(((BitSet) node.getMetaData("co")).get(i)){	
					    			throw new IllegalArgumentException("redo implementation");
//					    			co_seg = co_seg + "," + ancestralTreeLikelihoodInput.get().get(0).networkInput.get().segmentNames[i];
					    		}
					    	}
					    	buf.append(co_seg.replace("rem,", "") + "}");
					    	
					    	buf.append(",not={");
					    	String not_seg = "rem";
					    	for (int i = 0; i < ((BitSet) node.getMetaData("co")).size(); i++) {
					    		if(((BitSet) node.getMetaData("not")).get(i)){
					    			throw new IllegalArgumentException("redo implementation");
//					    			not_seg = not_seg + "," + ancestralTreeLikelihoodInput.get().get(0).networkInput.get().segmentNames[i];
					    		}
					    	}					    	
					    	buf.append(not_seg.replace("rem,", "") + "}");
					    	// check how "far" away from the exctinct corresponding lineage this one is
					    	for (Node l : ancestralTreeLikelihoodInput.get().get(0).tree.getExternalNodes()) {
					    		if (l.getID().startsWith(node.getID().replace("Survived", "Extinct"))) {
					    			int[] otherStates = new int[currSiteStates.length];
									for (int i = 0; i < ancestralTreeLikelihoodInput.get().size(); i++) {
									    int [] patternstates = ancestralTreeLikelihoodInput.get().get(i).getStatesForNode(
									    		ancestralTreeLikelihoodInput.get().get(0).tree, l);
									    
									    for (int j = 0; j < mapping.get(i).length; j++) {
									    	otherStates[mapping.get(i)[j]] = patternstates[j];
									    }	    
									}
									// compute the diffs
									int reassortmentDifference = 0;
									for (int i = 0; i < otherStates.length; i++) {
										if (otherStates[i]!=currSiteStates[i])
											reassortmentDifference++;
									}
							    	buf.append(",readist=" + reassortmentDifference);

					    		}
					    		
					    	}
					    	buf.append(",height=" + node.getHeight());


					    	
					    	
				    	}else if (node.getID().startsWith("reaExtinct_")) {
				    		// compute distance to the last node that has a distance measure, i.e. the last node that wasn't purely informed by the prior
				    		double dist = getLastNonPrior(node);	
					    	buf.append(",priordist=" + dist);
				    	}
			    	}
				    buf.append("]");
	    		}
	    	}
	    	
	    } else {
		    String seq_current = dataType.encodingToString(currSiteStates);
	    	if (logFullSequenceInput.get())  {
		    	buf.append("[&" + tagInput.get() + "=\"" + seq_current + "\"]");
	    	}else {
	    		if (!node.isRoot()) {
				    String seq_parent = dataType.encodingToString(parentSiteStates);
				    String muts = "remove";
				    for (int i = 0; i < seq_parent.length(); i++) {
				    	if (seq_parent.charAt(i) != seq_current.charAt(i)) {
				    		muts = muts + ","+seq_parent.charAt(i)+ ""+(i+1)+""+seq_current.charAt(i);
				    	}
				    }
				    muts = muts.replace("remove,", "");
				    muts = muts.replace("remove", "");
				    if (forTreeioInput.get())
				    	buf.append("[&" + tagInput.get() + "=\"" + muts.replace(",", "|") + "\"]");
				    else
				    	buf.append("[&" + tagInput.get() + "=\"" + muts + "\"]");

	    		}else {
			    	buf.append("[&" + tagInput.get() + "=\"" + seq_current + "\"]");

	    		}
	    	}
	    }

	    buf.append(':');
        buf.append(Math.max(node.getLength(),0.000000000001));
        return buf.toString();
    }
    
	private double getLastNonPrior(Node node) {
		if (node.getParent().getMetaData("distance")!=null)
			return node.getParent().getHeight()-node.getHeight();
		else
			return node.getParent().getHeight()-node.getHeight() + getLastNonPrior(node.getParent());
	}

	private void translate(int[] currSiteStates, List<List<String>> translated) {
		for (int a = 0; a < from.length; a++) {
			List<String> translatedFrame = new ArrayList<>();
			for (int i = from[a]; i < to[a]-2; i=i+3) {
				String codon = dataType.getCharacter(currSiteStates[i]) +""+dataType.getCharacter(currSiteStates[i+1])+""+dataType.getCharacter(currSiteStates[i+2]);
				translatedFrame.add(
						translationsTable.get(codon));
			}
			translated.add(translatedFrame);
		}
	}


	@Override
	public void close(PrintStream out) {
		((Tree) ancestralTreeLikelihoodInput.get().get(0).tree).close(out);
	}

	@Override
	public int getDimension() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double getArrayValue(int dim) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	private void getTOS(Tree tree){
		
		for (Node l : tree.getExternalNodes())
			if (!l.getID().startsWith("rea"))
				labelDist(l, 0.0);	
	}
	
	private void labelDist(Node n, double dist){
		if (n.isRoot())
			return;
		
		dist += n.getParent().getHeight()-n.getHeight();
		
		if(n.getParent().getMetaData("distance")==null)
			n.getParent().setMetaData("distance", dist);
		else
			n.getParent().setMetaData("distance", Math.max(dist, (double) n.getParent().getMetaData("distance")));
		
		// check if one of the child nodes has a 0 length
		if (n.getParent().getHeight() - n.getParent().getChild(0).getHeight()==0)
			n.getParent().getChild(0).setMetaData("distance", n.getParent().getMetaData("distance"));
		// check if one of the child nodes has a 0 length
		if (n.getParent().getHeight() - n.getParent().getChild(1).getHeight()==0)
			n.getParent().getChild(1).setMetaData("distance", n.getParent().getMetaData("distance"));

		labelDist(n.getParent(), dist);
	}


	
	
	private void initTable() {
		String dna = "AAA AAC AAG AAT ACA ACC ACG ACT AGA AGC AGG AGT ATA ATC ATG ATT CAA CAC CAG CAT CCA CCC CCG CCT CGA CGC CGG CGT CTA CTC CTG CTT GAA GAC GAG GAT GCA GCC GCG GCT GGA GGC GGG GGT GTA GTC GTG GTT TAA TAC TAG TAT TCA TCC TCG TCT TGA TGC TGG TGT TTA TTC TTG TTT";
		String aa = "K N K N T T T T R S R S I I M I Q H Q H P P P P R R R R L L L L E D E D A A A A G G G G V V V V O Y O Y S S S S O C W C L F L F";		
		
		translationsTable = new HashMap<>();
		String[] dnaarray = dna.split(" ");
		String[] aaarray = aa.split(" ");
		for (int i =0; i < dnaarray.length;i++)
			translationsTable.put(dnaarray[i], aaarray[i]);		
		
		
	}
}
