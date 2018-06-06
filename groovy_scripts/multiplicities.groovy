/*
 * author Timothy B. Hayward
 * May 2018
 * SIDIS dihadron (dipion) multiplicities
 * extract multiplicities from SIDIS CLAS12 events
 */

import java.io.*;

import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

import javax.swing.JFrame;
import org.jlab.groot.data.*;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.io.base.DataEvent;

import org.jlab.groot.math.F1D;
import org.jlab.groot.math.Func1D;
import org.jlab.groot.fitter.DataFitter;

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*; 
import histogram_calculations.*;

public class multiplicities {

	public static boolean dis_test(LorentzVector beam, LorentzVector target, 
			Particle pipi, Particle electron) {
    	double Q2 = 4*beam.e()*electron.e()*Math.pow(electron.theta()/2,2);

    	// target.e() is target mass in this case
    	double first_term = target.e()*target.e();
		double second_term = 2*target.e()*(beam.e()-electron.e());
		double third_term = - 4*beam.e()*electron.e()*Math.pow(Math.sin(electron.theta()/2),2);
		double W = Math.pow(first_term + second_term + third_term, 0.5);

		if ((Q2>1.0)&(W>2.0)) {
			return true;
		} else {
			return false;
		}
	}

	public static boolean sidis_test(LorentzVector beam, LorentzVector target, 
			Particle pipi, Particle electron) {
    	LorentzVector missing_part = new LorentzVector(pipi.px()+electron.px(), 
    		pipi.py()+electron.py(), beam.pz()-pipi.pz()-electron.pz(), 
    		beam.e()+target.e()-pipi.e()-electron.e());
    	double missing_mass = Math.pow(Math.pow(missing_part.e(),2)-
    		Math.pow(missing_part.p(),2),0.5);

		if (missing_mass>1.05) {
			return true;
		} else {
			return false;
		}
	}

	public static int array_index(def array, double value) {
		int index = 0;
		boolean index_exceeded = true;
		while(index_exceeded) {
			if (value < array[index]) {
				index_exceeded = false;
			} else {
				index++
			}
		}
		return (index-1);
	}

	private static double round (double value, int precision) {
		// just round a number to a certain precision
    	int scale = (int) Math.pow(10, precision);
    	return (double) Math.round(value * scale) / scale;
	}

	private static double mean (double x, double y) {
		return (x+y)/2;
	}

	private static double ratio_error (double ndh, ndis) {
		return Math.pow(ndh*(ndh+ndis)/Math.pow(ndis,3),0.5);
	}

	public static void main(String[] args) {
		// File[] hipo_list;
		// if (args.length == 0) {
		// 	// exits program if input directory not specified 
  //   	   	println("ERROR: Please enter a hipo file directory as the first argument");
  //   	  	System.exit(0);
  //   	} else {
  //   		File directory = new File(args[0]);
  //   		hipo_list = directory.listFiles();
  //   	}

        File root = new File("/volatile/clas12/data/rg-a/pass0/tag5b.3.3/");
        FilenameFilter beginswithm = new FilenameFilter() {
        	public boolean accept(File directory, String filename) {
         		return filename.startsWith("out_clas_003222");
         	}
        }
       	File[] hipo_list = root.listFiles(beginswithm);

        String[] run_list = ["out_clas_003973", "out_clas_003975"];
        // File root = new File("/volatile/clas12/data/rg-a/pass0/tag5b.3.3/");
        for (int i=0; i<run_list.size(); i++) {
        	FilenameFilter beginswith = new FilenameFilter() {
        		public boolean accept(File directory, String filename) {
         			return filename.startsWith(run_list[i]);
         		}
        	}
        	hipo_list = hipo_list+ root.listFiles(beginswith);
        }

       	for (File f: hipo_list) {
           	System.out.println(f);
        }

    	int n_files;
		if ((args.length < 2)||(Integer.parseInt(args[1])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			n_files = hipo_list.size();
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[1]);
		}

		// bin edges taken from proposal 
		// https://www.jlab.org/exp_prog/proposals/14/E12-06-112B_E12-09-008B.pdf
		def bins_Q2 =  [1.0, 1.32, 1.65, 2.09, 2.86, 11.0];
		def bins_xb = [0.0, 0.08, 0.1, 0.12, 0.14, 0.16, 0.19, 
			0.22, 0.26, 0.31, 0.8];
		def bins_z = [0.0, 0.31, 0.37, 0.43, 0.48, 0.53, 0.58, 0.63, 
			0.69, 0.76, 1.00]
		def bins_mpipi = [0.0, 0.38, 0.44, 0.50, 0.56, 0.64, 0.72, 
			0.78, 0.86, 1.06, 3.00];
		// def bins_Q2 =  [1.0, 1.32, 2.00, 2.86, 11.0];
		// def bins_xb = [0.0, 0.08, 0.1, 0.14, 0.20, 0.31, 0.8];
		// def bins_z = [0.0, 0.37, 0.48, 0.58, 0.69, 0.76, 1.00]
		// def bins_mpipi = [0.0, 0.38, 0.50, 0.64, 0.98, 3.00];

		int[][] counts_dis = new int[bins_Q2.size()-1][bins_xb.size()-1];
		int[][][][] counts_sidis = new int[bins_Q2.size()-1][bins_xb.size()-1][
			bins_z.size()-1][bins_mpipi.size()-1];

		double beam_mass = 0.00051099894; // electron mass (GeV)
		double target_mass = 0.93827208; // proton mass (GeV)
		LorentzVector beam = new LorentzVector(0, 0, 
			Math.pow(Math.pow(10.6,2)-Math.pow(beam_mass,2),0.5), 10.6); // 10.6 GeV beam
		LorentzVector target = new LorentzVector(0, 0, 0, target_mass);

		// GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6);
		GenericKinematicFitter rec_fitter = new event_builder_fitter(10.6);
		GenericKinematicFitter research_fitter = new enhanced_pid_fitter(10.6);
		EventFilter dis_filter = new EventFilter("11:X+:X-:Xn");
		EventFilter sidis_filter = new EventFilter("11:211:-211:X+:X-:Xn"); 

		int event_counter = 0;
		int dis_counter = 0;
		int sidis_counter = 0;
		for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files);
			// limit to a certain number of files defined by n_files
			HipoDataSource reader = new HipoDataSource();
			reader.open(hipo_list[current_file]);

			while(reader.hasEvent()==true){ // cycle through events
				event_counter++;
				HipoDataEvent event = reader.getNextEvent(); 
    			// PhysicsEvent  mc_Event  = mc_fitter.getPhysicsEvent(event);
    			PhysicsEvent  rec_Event  = rec_fitter.getPhysicsEvent(event);
    			PhysicsEvent  research_Event  = research_fitter.getPhysicsEvent(event);

    			if(dis_filter.isValid(research_Event)==true){
    				Particle pipi = research_Event.getParticle("[211]+[-211]");
    				Particle electron = research_Event.getParticle("[11]");

    				if (dis_test(beam, target, pipi, electron)) {
    					dis_counter++;
    					double Q2 = 4*beam.e()*electron.e()*Math.pow(electron.theta()/2,2);
    					double xb = Q2 / (2*target.e()*(beam.e()-electron.e()));

    					if ((Q2<11.0)&(xb<0.8)) {
    						counts_dis[array_index(bins_Q2, Q2)][array_index(bins_xb, xb)]++;
    					}

    					if(sidis_filter.isValid(research_Event)==true){
    						if (sidis_test(beam, target, pipi, electron)) {
    							sidis_counter++;
    							double z = pipi.e()/(beam.e()-electron.e());
    							double mpipi = pipi.mass();
    							if (((Q2<11.0)&(xb<0.8))&((z<1.0)&(mpipi<3.00))) {
    								counts_sidis[array_index(bins_Q2, Q2)][array_index(
    									bins_xb, xb)][array_index(bins_z, z)][array_index(
    									bins_mpipi, mpipi)]++;
    							}


    						}
    					}
    				}
    			}
    		}
		}

		println(counts_dis);
		println(counts_sidis);
		println();
		println("Analyzed "+event_counter+" events.");
		println(dis_counter+" dis events ("+
			round(100*dis_counter/event_counter,3)+"%)");
		println(sidis_counter+" sidis events ("+
			round(100*sidis_counter/event_counter,3)+"%)");


		// JFrame[] frame = new JFrame[bins_z.size()-1];
		// EmbeddedCanvas canvas_mpipi = new EmbeddedCanvas();

		// for (int j=0; j < bins_z.size()-1; j++) {
		// 	frame[j].setSize(1200,900);
		// 	GraphErrors ratio_points_mpipi = new GraphErrors();
		// 	for (int i=0; i<bins_mpipi.size()-1; i++) {
		// 		if ((counts_dis[0][0]>0)&(counts_sidis[0][0][j][i]>0)) {
		// 			double bin_center = mean(bins_mpipi[i],bins_mpipi[i+1]);
		// 			double ratio = counts_sidis[0][0][j][i]/counts_dis[0][0];
		// 			double error = ratio_error(counts_sidis[0][0][j][i],counts_dis[0][0]);
		// 			ratio_points_mpipi.addPoint(bin_center, ratio, 0, error);
		// 		}
		// 	}
		// 	ratio_points_mpipi.setMarkerColor(0);
		// 	ratio_points_mpipi.setLineColor(0);
		// 	ratio_points_mpipi.setMarkerStyle(0);
		// 	ratio_points_mpipi.setMarkerSize(4);
		// 	ratio_points_mpipi.setLineThickness(2);
		// 	ratio_points_mpipi.setTitleX("M#pi#pi (GeV)");
		// 	ratio_points_mpipi.setTitleY("Multiplicity");
		// 	// canvas_mpipi.getPad(j).setTitle("Hadron Pair Fractional Virtual Photon Energy");
		// 	canvas_mpipi.setTitleSize(20);
		// 	canvas_mpipi.setAxisTitleSize(18);
		// 	canvas_mpipi.setAxisLabelSize(18);
		// 	canvas_mpipi.setStatBoxFontSize(18);
		// 	canvas_mpipi.draw(ratio_points_mpipi);
		// 	frame.add(canvas_mpipi);
		// 	frame.setLocationRelativeTo(null);
		// 	frame.setVisible(true);
		// }


		PrintStream file_counts_dis = new PrintStream(new File("counts_dis.txt"));
		// Store current System.out before assigning a new value
        PrintStream console = System.out;
        // Assign file_counts_dis to output stream
        System.setOut(file_counts_dis);
        for (int i = 0; i < bins_Q2.size()-1; i++) {
        	String line = "";
        	for (int j = 0; j < bins_xb.size()-1; j++) {
        		line+=" ";
        		line+=Integer.toString(counts_dis[i][j]);
        	}
        	System.out.println(line);
        }

        PrintStream file_counts_sidis = new PrintStream(new File("counts_sidis.txt"));
        // Assign file_counts_sidis to output stream
        System.setOut(file_counts_sidis);
        for (int i = 0; i < bins_Q2.size()-1; i++) {
        	for (int j = 0; j < bins_xb.size()-1; j++) {
        		for (int n = 0; n < bins_z.size()-1; n++) {
        			String line = ""
        			for (int m = 0; m < bins_mpipi.size()-1; m++) {
        				line+=" ";
        				line+=Integer.toString(counts_sidis[i][j][n][m]);
        			}
        			System.out.println(line);
        		}
        		// System.out.println();
        	}
        	// System.out.println();
        }

        PrintStream file_bins = new PrintStream(new File("bins.txt"));
        //Assign file_bins to output stream
        System.setOut(file_bins);
        String line = ""
        for (int i = 0; i < bins_Q2.size(); i++) {
        	line+=" "+Double.toString(bins_Q2[i]);
        }
        System.out.println(line);
        line = ""
        for (int i = 0; i < bins_xb.size(); i++) {
        	line+=" "+Double.toString(bins_xb[i]);
        }
        System.out.println(line);
        line = ""
        for (int i = 0; i < bins_z.size(); i++) {
        	line+=" "+Double.toString(bins_z[i]);
        }
        System.out.println(line);
        line = ""
        for (int i = 0; i < bins_mpipi.size(); i++) {
        	line+=" "+Double.toString(bins_mpipi[i]);
        }
        System.out.println(line);


	}
}