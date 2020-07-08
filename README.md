# LDscaff
LDscaff: LD-based scaffolding of de novo genome assemblies


We proposed LDscaff for draft genome assembly incorporating linkage disequilibrium information in population sequencing data. LDscaff aids to find a optimal layout of pre-assembled scaffolds with a graph method. By taking as input the population variation data and scaffolds, LDscaff builds an undirected graph with vertices and edges, representing the scaffold ends and the LD strength among them. Computing the optimal orders and orientations of these scaffolds can be treated as a maximum weight matching(MWM) problem. After solving the MWM problem, users can get a updated re-assembled scaffolds.
