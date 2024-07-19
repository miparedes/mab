package nab.util;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.Sequence;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.datatype.DataType;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.BEASTClassLoader;
import beast.util.PackageManager;
import beast.util.Randomizer;
import feast.nexus.CharactersBlock;
import feast.nexus.NexusBuilder;
import feast.nexus.TaxaBlock;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Tim Vaughan
 */
@Description("An adapted version of SimulatedAlignment for TGV to simulate sequences with dnds ratios.")
public class SimulatedDNDSAlignment extends Alignment {

    public Input<Tree> treeInput = new Input<>(
            "tree",
            "Tree down which to simulate sequence evolution.",
            Input.Validate.REQUIRED);

    public Input<SiteModel> siteModelInput = new Input<>(
            "siteModel",
            "Site model to use in simulation.",
            Input.Validate.REQUIRED);

    public Input<Integer> sequenceLengthInput = new Input<>(
            "sequenceLength",
            "Length of sequence to simulate.",
            Input.Validate.REQUIRED);

    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName",
            "Name of file (if any) simulated alignment should be saved to.");

    private Tree tree;
    private SiteModel siteModel;
    private int seqLength;
    private DataType dataType;

    private String ancestralSeqStr;

    public SimulatedDNDSAlignment() {
        sequenceInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {

        tree = treeInput.get();
        siteModel = siteModelInput.get();
        seqLength = sequenceLengthInput.get();

        sequences.clear();

        grabDataType();

        simulate();

        super.initAndValidate();

        // Write simulated alignment to disk if required
        if (outputFileNameInput.get() != null) {
            try (PrintStream pstream = new PrintStream(outputFileNameInput.get())) {
                NexusBuilder nb = new NexusBuilder();
                nb.append(new TaxaBlock(new TaxonSet(this)));
                nb.append(new CharactersBlock(this));
                nb.write(pstream);
            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error writing to file "
                        + outputFileNameInput.get() + ".");
            }
        }
    }

    /**
     * Perform actual sequence simulation.
     */
    private void simulate() {
        int nTaxa = tree.getLeafNodeCount();

        double[] categoryProbs = siteModel.getCategoryProportions(tree.getRoot());

        int nCategories = siteModel.getCategoryCount();
        int nStates = dataType.getStateCount();
        double[][] transitionProbs = new double[nCategories][nStates*nStates];

        int[][] alignment = new int[nTaxa][seqLength];

        int[] categories = new int[seqLength];
        for (int i=0; i<seqLength; i++)
            categories[i] = Randomizer.randomChoicePDF(categoryProbs);

        Node root = tree.getRoot();

        int[] parentSequence = new int[seqLength];
        double[] frequencies = siteModel.getSubstitutionModel().getFrequencies();
        for (int i=0; i<parentSequence.length; i++)
            parentSequence[i] = Randomizer.randomChoicePDF(frequencies);

        ancestralSeqStr = dataType.encodingToString(parentSequence);

        traverse(root, parentSequence,
                categories, transitionProbs,
                alignment);

        for (int leafIdx=0; leafIdx<nTaxa; leafIdx++) {
            String seqString = dataType.encodingToString(alignment[leafIdx]);

            String taxonName;
            if (tree.getNode(leafIdx).getID() != null)
                taxonName = tree.getNode(leafIdx).getID();
            else
                taxonName = "t" + leafIdx;

            sequenceInput.setValue(new Sequence(taxonName, seqString), this);
        }
    }

    /**
     * Traverse a tree, simulating a sequence alignment down it.
     *
     * @param node Node of the tree
     * @param parentSequence Sequence at the parent node in the tree
     * @param categories Mapping from sites to categories
     * @param transitionProbs transition probabilities
     * @param regionAlignment alignment for particular region
     */
    private void traverse(Node node,
            int[] parentSequence,
            int[] categories, double[][] transitionProbs,
            int[][] regionAlignment) {

        for (Node child : node.getChildren()) {

            // Calculate transition probabilities
            for (int i=0; i<siteModel.getCategoryCount(); i++) {
                siteModel.getSubstitutionModel().getTransitionProbabilities(
                        child, node.getHeight(), child.getHeight(),
                        siteModel.getRateForCategory(i, child),
                        transitionProbs[i]);
            }

            // Draw characters on child sequence
            int[] childSequence = new int[parentSequence.length];
            int nStates = dataType.getStateCount();
            double[] charProb = new double[nStates];
            for (int i=0; i<childSequence.length; i++) {
                int category = categories[i];
                System.arraycopy(transitionProbs[category],
                        parentSequence[i]*nStates, charProb, 0, nStates);
                childSequence[i] = Randomizer.randomChoicePDF(charProb);
            }

            if (child.isLeaf()) {
                System.arraycopy(childSequence, 0,
                        regionAlignment[child.getNr()], 0, childSequence.length);
            } else {
                traverse(child, childSequence,
                        categories, transitionProbs,
                        regionAlignment);
            }
        }
    }

    /**
     * HORRIBLE function to identify data type from given description.
     */
    private void grabDataType() {
        if (userDataTypeInput.get() != null) {
            dataType = userDataTypeInput.get();
        } else {

            List<String> dataTypeDescList = new ArrayList<>();
            List<String> classNames = PackageManager.find(beast.evolution.datatype.DataType.class, "beast.evolution.datatype");
            for (String className : classNames) {
                try {
                    DataType thisDataType = (DataType) BEASTClassLoader.forName(className).newInstance();
                    if (dataTypeInput.get().equals(thisDataType.getTypeDescription())) {
                        dataType = thisDataType;
                        break;
                    }
                    dataTypeDescList.add(thisDataType.getTypeDescription());
                } catch (ClassNotFoundException
                    | InstantiationException
                    | IllegalAccessException e) {
                }
            }
            if (dataType == null) {
                throw new IllegalArgumentException("Data type + '"
                        + dataTypeInput.get()
                        + "' cannot be found.  Choose one of "
                        + Arrays.toString(dataTypeDescList.toArray(new String[0])));
            }
        }
    }

}
