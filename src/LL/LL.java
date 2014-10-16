/*
 * Name: Cedirc Jo
 * Project 6
 * 
 */

package LL;

import java.io.*;
import java.util.*;
import java.text.*;


public class LL {

    private int column;
    private int row;
    private double[][] data;
    private int k;  // number of clusters
    ArrayList<Integer> r;   // set of random numbers for generating cluster center
    private int rnd;    //
    private int[][] cluster;
    private int loop;   // number of iteration
    
    ArrayList<Double> cluster1, cluster2, cluster3;
    int count1, count2, count3;
    
    
    public LL (int column, int row, double[][] data, int k) {
        this.column = column;
        this.row = row;
        this.data = data;
        this.k = k;
        this.cluster = new int[row][column];
        this.loop = 1;
        this.cluster1 = new ArrayList<>();
        this.cluster2 = new ArrayList<>();
        this.cluster3 = new ArrayList<>();
        count1 = 0;
        count2 = 0;
        count3 = 0;
       
        
        this.r = new ArrayList<>();
        for(int i=1; i <= row; i++) {
            r.add(i);
        }
        Collections.shuffle(r);
        
        rnd = (int)(Math.random() * 11);
        
        lloyd(cluster1, cluster2, cluster3, cluster);
    }
    
    
    public void lloyd(ArrayList cst1, ArrayList cst2, ArrayList cst3, int[][] cluster) {
        int[][] clust = cluster;
        /*--------------------------------------------------------------------*/
        // Pick random points
        if (k >= 1) {
            int t = random();
            for(int i=0; i < column; i++) {
                cst1.add(data[t][i]);
            }
        }
         
        if (k >= 2) {
            int t = random();
            for(int i=0; i < column; i++) {
                cst2.add(data[t][i]);
            }
        }
        
        if (k >= 3) {
            int t = random();
            for(int i=0; i < column; i++) {
                cst3.add(data[t][i]);
            }
        }
        
        /*--------------------------------------------------------------------*/
        while (!terminate(clust, setCluster(cst1, cst2, cst3))) {
            clust = setCluster(cst1, cst2, cst3);
         
            /*--------------------------------------------------------------------*/
            // Update to new cluster center
            cst1 = newCenter2(cst1, count1, 1, clust);
            cst2 = newCenter2(cst2, count2, 2, clust);
            cst3 = newCenter2(cst3, count3, 3, clust);
                   
            setCluster(cst1, cst2, cst3);
            
            loop++;
        }
        
/*----------------------------------------------------------------------------*/
        // Output
        System.out.print("Genes in Cluster 1: ");
        for (int i=0; i < count1; i++) {
            int t = clust[i][0];
            System.out.print("g" + t + " ");
        }
        System.out.println();

        System.out.print("Genes in Cluster 2: ");

        for (int i=0; i < count2; i++) {
            int t = clust[i][1];
            System.out.print("g" + t + " ");
        }
        System.out.println();

        System.out.print("Genes in Cluster 3: ");
        for (int i=0; i < count3; i++) {
            int t = clust[i][2];
            System.out.print("g" + t + " ");
        }
        System.out.println();

        System.out.print("Cluster 1 Center: ");
        for (int j=0; j < k; j++) {
            System.out.print(cst1.get(j) + " ");
        }
        System.out.println();

        System.out.print("Cluster 2 Center: ");
        for (int j=0; j < k; j++) {
            System.out.print(cst2.get(j) + " ");
        }
        System.out.println();

        System.out.print("Cluster 3 Center: ");
        for (int j=0; j < k; j++) {
            System.out.print(cst3.get(j) + " ");
        }
        System.out.println();
        

        double sed = sed(cst1, clust, count1, 1) + 
                     sed(cst2, clust, count2, 2) +
                     sed(cst3, clust, count3, 3);
        DecimalFormat df = new DecimalFormat("#.##");
        System.out.println("Squared error distortion: " + df.format(sed/k));

        System.out.println("Number of iterations: " + loop);  
    }
    
/******************************************************************************/
    //squared error distortion
    public double sed(ArrayList<Double> cst, int[][] clust, int count, int k) {
        double result = 0;
        double tmp = 0;
        
        for(int i=0; i < count; i++) {
            int loc = clust[i][k-1]-1;
  
            
            for(int j=0; j < column; j++) {
                tmp = cst.get(j) - data[loc][j];
                result += Math.pow(tmp,2);
            }
        }
        return result;
    }
    
/******************************************************************************/
    //Generate a random point
    public int random() {
        rnd = (rnd+1) % r.size();
        //System.out.println("rnd = " + rnd);
        return rnd; 
    }
    
    
/******************************************************************************/
    //Returns Euclidean distance between two points
    // a = cluster point
    public double distance(ArrayList<Double> a, int b) {
        double tmp1, tmp2;
        tmp2 = 0;
        
        for(int i=0; i < column; i++) {
            tmp1 = a.get(i) - data[b][i];
            tmp2 += Math.pow(tmp1, 2);
        }
        return Math.sqrt(tmp2);
    }
    
    
/******************************************************************************/    
    //Returns a list of all distances between random point and all data points
    public ArrayList distBtwPts(ArrayList a) {
        ArrayList<Double> dist = new ArrayList<>();
        
        for (int i=0; i < row; i++) {
            double d = distance(a, i);
            dist.add(d);
        }
        return dist;   
    }
    
/******************************************************************************/    
    //Update to new cluster center    
    public ArrayList newCenter2(ArrayList c, int count, int k, int[][] clust) {
        ArrayList<Double> cc = new ArrayList<>();
        
        double xMean = 0, yMean = 0, zMean = 0;
        
        for (int i=0; i < count; i++) {
            if (clust[i][k-1]-1 > 0) { 
            int x = clust[i][k-1]-1;
             xMean += data[x][0];
             yMean += data[x][1];
             zMean += data[x][2];
        }}
        DecimalFormat df = new DecimalFormat("#.#");
        xMean = Double.parseDouble(df.format(xMean/count));
        yMean = Double.parseDouble(df.format(yMean/count));
        zMean = Double.parseDouble(df.format(zMean/count));
        
        cc.add(xMean);
        cc.add(yMean);
        cc.add(zMean);
         
        return cc;
    }
    
    
    
    
/******************************************************************************/    
    //Group genes into each cluster
    public int[][] setCluster(ArrayList cst1, ArrayList cst2, ArrayList cst3) {
        int[][] sc = new int[row][column];
        this.count1 = 0;    
        this.count2 = 0;
        this.count3 = 0;
        
        ArrayList<Double> c1 = distBtwPts(cst1);
        ArrayList<Double> c2 = distBtwPts(cst2);
        ArrayList<Double> c3 = distBtwPts(cst3);
       
        
        for (int j=0; j < row; j++) {
            double min = Math.min(c3.get(j), Math.min(c1.get(j), c2.get(j)));
            
            if(min == c1.get(j)) {
                sc[count1][0] = j+1;
                count1++;
            }
            else if(min == c2.get(j)) {
                sc[count2][1] = j+1;
                count2++;
            }
            else if(min == c3.get(j)) {
                sc[count3][2] = j+1;
                count3++;
            }
        }
        return sc;
    }
    
    
  
    public boolean terminate (int[][] a, int[][] b) {
        boolean tf = true;

        for(int i=0; i < row; i++) {
            for(int j=0; j < column; j++) {
                if (a[i][j] != b[i][j]) {
                    tf = false;
                    break;
                }
            }
        }
        return tf;
    }
    
    

/*============================================================================*/
    
    
public static void main (String[] args) throws IOException {
    
    int row = 0;
    int column = 0;
    String filename = "im.txt";
    Scanner input = new Scanner(new File(filename));
    
    /*------------------------------------------------------------------------*/
    /*  Read the number of genes and the number of time points
     */
    row = Integer.parseInt(input.next());
    column = Integer.parseInt(input.next());
    //System.out.print("column: " + column + " row: " + row);
    /*------------------------------------------------------------------------*/
    
    ArrayList<Double> list = new ArrayList<>();
    
    /*------------------------------------------------------------------------*/
    /*  Read the intensity at each time point
     */
    while(input.hasNext())
    {
        String line = input.next();
        
        for (int i=0; i < column; i++) {
            String in = input.next();
            list.add(Double.parseDouble(in));
        }
    }
    
    /*------------------------------------------------------------------------*/
    double[][] data = new double[row][column];
    
    int k=0;
    while (k < list.size()) {
        for(int i=0; i < row; i++) {
            for(int j=0; j < column; j++) {
                data[i][j] = list.get(k);
                k++;
            }
        }
    }
    
    
//    System.out.print("Enter the number of clusters: ");
//    String in = new Scanner(System.in).next();
//    int c = Integer.parseInt(in);
    int c = 3;
    
    LL ll = new LL(column, row, data, c);
    }
}
