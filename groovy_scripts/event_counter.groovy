/*
 * author Timothy B. Hayward
 * April 2018
 * SIDIS dihadron (dipion) event counter wtih various fitters
 */

import java.io.File;

import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

import javax.swing.JFrame;
import org.jlab.groot.data.*;
import org.jlab.groot.graphics.EmbeddedCanvas;
import org.jlab.io.base.DataEvent;

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*; // 

public class event_counter {
	// groovy program to calculate the rate of reconstructed dipion events following CLAS12 
	// detector simulation and reconstruction (gemc / myClara)
	// c.f. /volatile/clas12/thayward/SIDIS/dihadron/README

	private static double round (double value, int precision) {
		// just round a number to a certain precision
    int scale = (int) Math.pow(10, precision);
    return (double) Math.round(value * scale) / scale;
	}

	public static void main(String[] args) {

		File[] hipo_list;
		if (args.length == 0) {
			// exits program if input directory not specified 
        	println("ERROR: Please enter a hipo file directory as the first argument");
       		System.exit(0);
    	} else {
    		File directory = new File(args[0]);
    		hipo_list = directory.listFiles();
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

		GenericKinematicFitter mc_fitter = new monte_carlo_fitter(10.6);
		GenericKinematicFitter generic_fitter = new generic_rec_fitter(10.6);
		// GenericKinematicFitter research_fitter = new SIDIS_fitter(10.6);
		EventFilter filter = new EventFilter("-211:X+:X-:Xn"); 
			// all events with electron, + and - pion and any number of other particles

		int event_counter = 0; // number of events analyzed
		int monte_carlo_counter = 0; // number of events with reconstructed MC particles
		int generic_counter = 0; // number of events reconstructed from generic CLAS12 fitter
		// int research_counter = 0; // number of events reconstructed with research project fitter

		for (int current_file; current_file<n_files; current_file++) {
		// for (int current_file; current_file<n_files; current_file++) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files);
			// limit to a certain number of files defined by n_files
			HipoDataSource reader = new HipoDataSource();
			reader.open(hipo_list[current_file]); // open next hipo file

			while(reader.hasEvent()==true){ // cycle through events
				event_counter++; // new event analyzed
				HipoDataEvent event = reader.getNextEvent(); 
    			PhysicsEvent  mc_Event  = mc_fitter.getPhysicsEvent(event);
    			PhysicsEvent  generic_Event  = generic_fitter.getPhysicsEvent(event);
    			// PhysicsEvent  research_Event  = research_fitter.getPhysicsEvent(event);

    			if(filter.isValid(mc_Event)==true){
    				monte_carlo_counter++; // monte carlo fitter returned recon dipion event
    			}
    			if(filter.isValid(generic_Event)==true){
    				generic_counter++; // generic CLAS12 fitter returned recon dipion event
    			}
    			// if(filter.isValid(research_Event)==true){
    			// 	research_counter++; // research fitter returned recon dipion event
    			// }
			}
		}
		println(); println();
		println("Analyzed "+event_counter+" events.");
		println("monte_carlo_fitter reconstructed "+monte_carlo_counter+" events ("+
			round(100*monte_carlo_counter/event_counter,3)+"%)");
		println("generic fitter reconstructed "+generic_counter+" events ("+
			round(100*generic_counter/event_counter,3)+"%)");
		// println("research project fitter reconstructed "+research_counter+" events ("+
		// 	round(100*research_counter/event_counter,3)+"%)");
		println();
	}
}
