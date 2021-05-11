/*
 * author Timothy B. Hayward
 * 
 * SIDIS dihadron 
 */

import java.io.File;

import org.jlab.io.hipo.*;
import org.jlab.io.base.DataEvent;
import org.jlab.clas.physics.*;
import org.jlab.clas12.physics.*;

// import from hayward_coatjava_extensions
import extended_kinematic_fitters.*; 
import analyzers.*;

// dilks CLAS QA analysis
import clasqa.QADB

// filetype for gathering files in directory
import groovy.io.FileType;


public class processing_dihdarons {

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

    	def list = []

		def dir = new File(args[0])
		dir.eachFileRecurse (FileType.FILES) { file ->
		  list << file
		  // println(file.toString()); println(); 
		}

		println(); println(); println();

		String p1_Str;
		if (args.length < 2) {
			// assigns pi+ to p1
			println("WARNING: Specify a PDG PID for p1! Set to pi+. \n");
			p1_Str = "211";
		} else {
			p1_Str = args[1];
			println("Set p1 PID = "+p1_Str+"\n");
		}

		String p2_Str;
		if (args.length < 3) {
			// assigns pi- to p2
			println("WARNING: Specify a PDG PID for p2! Set to pi-. \n");
			p2_Str = "-211";
		} else {
			p2_Str = args[2];
			println("Set p2 PID = "+p2_Str+"\n");
		}

		String output_file;
		if (args.length < 4) {
			// uses dummy name for output file if not specified
			println("WARNING: Specify an output file name. Set to \"dihadron_dummy_out.txt\".\n");
			output_file = "dihadron_dummy_out.txt"
		} else {
			output_file = args[3];
		}

		int n_files;
		if ((args.length < 5)||(Integer.parseInt(args[4])>hipo_list.size())) {
			// if number of files not specified or too large, set to number of files in directory
			println("WARNING: Number of files not specified or number too large."); 
			println("Setting # of files to be equal to number of files in the directory.");
			// n_files = hipo_list.size();
			n_files = list.size();
			println("There are "+hipo_list.size()+" or maybe "+list.size()+" number of files.")
		} else{
			// if specified, convert to int
			n_files = Integer.parseInt(args[4]);
		}

		File file = new File(output_file);
		file.bytes = new byte[0]

		int hadron_pair_counts = 0;
		GenericKinematicFitter research_fitter = new analysis_fitter(10.6041); // load my kinematic fitter/PID
		EventFilter filter = new EventFilter("11:"+p1_Str+":"+p2_Str+":X+:X-:Xn"); // set filter for final states
		// setup QA database
		QADB qa = new QADB();

		int num_events = 0;
		int current_file = 0;
		// for (int current_file; current_file<n_files; current_file++) {
		while (current_file < n_files) {
			println(); println(); println("Opening file "+Integer.toString(current_file+1)
				+" out of "+n_files); println(); println();
			// limit to a certain number of files defined by n_files

			HipoDataSource reader = new HipoDataSource();

			reader.open(list[current_file]); // open next hipo file
			current_file++;
			HipoDataEvent event = reader.getNextEvent(); 

			while(reader.hasEvent()==true){
				num_events++; 
				if (num_events%250000 == 0) { // not necessary
					print("processed: "+num_events+" events. ");
				}

				// get run and event numbers
				event = reader.getNextEvent();
			    int runnum = event.getBank("RUN::config").getInt('run',0); // collect info for QA
			    int evnum = event.getBank("RUN::config").getInt('event',0);

			    PhysicsEvent research_Event = research_fitter.getPhysicsEvent(event);

			    boolean process_event = false;
			    if (runnum == 11) { // if run number = 11 then it is MC and we don't use QA
			    	process_event = filter.isValid(research_Event);
			    } else {
			    	process_event = (filter.isValid(research_Event) && qa.OkForAsymmetry(runnum,evnum));
			    }
				if (process_event) {
					int num_p1 = research_Event.countByPid(p1_Str.toInteger());  // get # of particles w/ pid1
					int num_p2 = research_Event.countByPid(p2_Str.toInteger()); // get # of particles w/ pid2

					for (int current_p1 = 0; current_p1 < num_p1; current_p1++) { // cycle over all combinations
						for (int current_p2 = 0; current_p2 < num_p2; current_p2++) {

							Dihadrons variables = new Dihadrons(event, research_Event, 
								p1_Str.toInteger(), current_p1, p2_Str.toInteger(), current_p2);
							// this is my class for defining all relevant kinematic variables

							if (variables.channel_test(variables)) {
								int helicity = variables.get_helicity(); // helicity of event, might be 0

								// lab kinematics
								double e_p = variables.e_p();
								double e_theta = variables.e_theta();
								double p1_p = variables.p1_p();
								double p1_theta = variables.p1_theta();
								double p2_p = variables.p2_p();
								double p2_theta = variables.p2_theta();

								// DIS variables
								double Q2 = variables.Q2();
								double W = variables.W();
								double y = variables.y();
								double Mx = variables.Mx();
								double Mx1 = variables.Mx1();
								double Mx2 = variables.Mx2();

								// SIDIS variables
								double x = variables.x();
								double z = variables.z();
								double xF = variables.xF();
								double pT = variables.pT();
								double eta = variables.eta();
								double eta_gN = variables.eta_gN();

								// SIDIS dihadron variables
								double z1 = variables.z1();
								double z2 = variables.z2();
								double xF1 = variables.xF1();
								double xF2 = variables.xF2();
								double zeta = variables.zeta(); // assumption here is p1 is a proton
								// otherwise zeta is meaningless
								double Mh = variables.Mh();
								double pT1 = variables.pT1();
								double pT2 = variables.pT2();
								double pTpT = variables.pTpT();
								double eta1 = variables.eta1();
								double eta2 = variables.eta2();
								double Delta_eta = variables.Delta_eta();
								double eta1_gN = variables.eta1_gN();
								double eta2_gN = variables.eta2_gN();

								// angles 
								double phi1 = variables.phi1();
								double phi2 = variables.phi2();
								double Delta_phi = variables.Delta_phi();
								double phih = variables.phih();
								double phiR = variables.phiR();
								double theta = variables.theta();

								// depolarization factors
								double Depolarization_A = variables.Depolarization_A();
								double Depolarization_B = variables.Depolarization_B();
								double Depolarization_C = variables.Depolarization_C();
								double Depolarization_V = variables.Depolarization_V();
								double Depolarization_W = variables.Depolarization_W();

								// append event to next line of the text file
								file.append(runnum+" "+evnum+" "+helicity+" ");
								file.append(e_p+" "+e_theta+" "+p1_p+" "+p1_theta+" "+p2_p+" "+p2_theta+" ");
								file.append(Q2+" "+W+" "+Mx+" "+Mx1+" "+Mx2+" "+x+" "+y+" "+z+" ");
								file.append(z1+" "+z2+" "+Mh+" "+xF+" "+xF1+" "+xF2+" ");
								file.append(pT+" "+pT1+" "+pT2+" "+pTpT+" "+zeta+" ");
								file.append(eta+" "+eta1+" "+eta2+" "+Delta_eta+" "+eta1_gN+" "+eta2_gN+" ");
								file.append(phi1+" "+phi2+" "+Delta_phi+" "+phih+" "+phiR+" "+theta+" ");
								file.append(Depolarization_A+" "+Depolarization_B+" "+Depolarization_C+" ");
								file.append(Depolarization_V+" "+Depolarization_W+"\n");
							}
						}
					}
				}
			}
			println(); println();
			print("1:runnum, 2:evnum, 3:helicity, 4:e_p, 5:e_theta, 6:p1_p, 7:p1_theta, ");
			print("8:p2_p, 9:p2_theta, 10:Q2, 11:W, 12:Mx, 13:Mx1, 14:Mx2, 15:x, 16:y, 17:z ");
			print("18:z1, 19:z2, 20:Mh, 21:xF, 22:xF1, 23:xF2, ");
			print("24:pT, 25:pT1, 26:pT2, 27:pTpT 28:zeta ");
			print("29:eta, 30:eta1, 31:eta2, 32:Delta_eta, 33:eta1_gN, 34:eta2_gN ");
			print("35:phi1, 36:phi2, 37:Delta_phi, 38:phih, 39:phiR, 40:theta");
			print("41:DepA, 42:DepB, 43:DepC, 44:DepV, 45:DepW, 46: MxRad, 47: Mx1Rad, 48: Mx2Rad. \n");

			println(); println();
			println("Set p1 PID = "+p1_Str+"\n");
			println("Set p2 PID = "+p2_Str+"\n");
			println("output file is: "+file);
		}

	}
}