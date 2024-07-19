/*
* File ClusterTree.java
*
* Copyright (C) 2010 Remco Bouckaert remco@cs.auckland.ac.nz
*
* This file is part of BEAST2.
* See the NOTICE file distributed with this work for additional
* information regarding copyright ownership and licensing.
*
* BEAST is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
*  BEAST is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public
* License along with BEAST; if not, write to the
* Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
* Boston, MA  02110-1301  USA
*/
package nab.util;


import beast.core.Description;
import beast.core.Input;
import beast.core.StateNodeInitialiser;
import beast.evolution.alignment.distance.Distance;
import beast.evolution.tree.Node;
import beast.util.ClusterTree;



/**
 * Adapted from Weka's HierarchicalClustering class *
 */
@Description("Resolves polytomies of cluster trees for initial states")
public class ResolvedClusterTree extends ClusterTree implements StateNodeInitialiser {

    final public Input<Double> heightInput = new Input<>("rootHeight", "scales the nodes to match the root height");

	
    @Override
    public void initAndValidate() {
    	super.initAndValidate();
    	double scaleFactor = heightInput.get()/getRoot().getHeight();
        for (Node node : getInternalNodes()) {
            double height = node.getHeight();
            node.setHeight(height*scaleFactor);
        }
    	System.out.println(toString());
//    	System.exit(0);
        initStateNodes();
    }

} // class ClusterTree
