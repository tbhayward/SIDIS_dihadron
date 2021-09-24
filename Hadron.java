/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package analyzers;

/**
 *
 * @author tbhayward
 */

import org.jlab.clas.physics.Particle;
import org.jlab.clas.physics.PhysicsEvent;
import org.jlab.io.base.DataEvent;
import org.jlab.io.hipo.HipoDataBank;
import org.jlab.clas.physics.*;


public class Hadron {
    
    protected byte helicity;
    protected int runnum;
    
    protected double test;
    
    protected int num_elec, num_piplus, num_piminus, num_kplus, num_kminus, num_protons, num_particles;
    
    // labels are unnumbered if they refer to the dihadron (perhaps a meson) and numbered for individual
    // hadrons. Convention is ordered by mass, then charge. For example in pi+pi- pi+ is hadron 1
    // in proton+pi+ the proton is p1, in k+pi- the kaon is p1.
    protected double Q2, W, gamma, nu, x, y, z;
    protected double Mx, Mx2; 
    protected double Mh, pT, xF, zeta;
    protected double eta, eta_gN;
    // eta is the rapidity, preferred by theorists in the Breit frame (e.g. eta1 is in Breit) 
    // eta_gN is the rapidity in the gamma*-nucleon COM frame
    // the difference between two rapidities is Lorentz invariant, i.e.
    // eta1-eta2 = eta1_COM - eta2_COM

    protected double phi;
    
    // depolarization vectors defining the polarization lost during the transfer from beam to 
    // the virtual photon. 
    // in ALU BSAs the twist 2 terms are modified by C/A and the twist 3 terms by W/A
    // B and V come in AUL
    protected double Depolarization_A;
    protected double Depolarization_B;
    protected double Depolarization_C;
    protected double Depolarization_V;
    protected double Depolarization_W;
    
    protected double e_px, e_py, e_pz, e_p, e_e, e_theta, e_phi; // electron kinematics
    protected double p_px, p_py, p_pz, p_p, p_e, p_theta, p_phi; // dihadron kinematics, mass is Mh
    protected double vz_e, vz_p;
    
    protected double p_Breit_pz, p_gN_pz;
    
    public static boolean channel_test(Hadron variables) {
//        if (variables.helicity==0){ 
//            System.out.println("You're returning false because helicity = 0. Is this data or MC?");
//            return false; }
        if (variables.Q2()<1) { return false; } 
        if (variables.W()<2) { return false; } 
//        else if (variables.xF1()<0.0 || variables.xF2()<0.0) { return false; } 
//        else if (variables.y()>0.80) { return false; } 
        else if (variables.Mx()<0.0) { return false; } 
	return true;
    }
    
    public Hadron(DataEvent event, PhysicsEvent recEvent, int pPID, int pIndex) {
        // provide the PDG PID of the two hadrons
        
        // load banks
        HipoDataBank eventBank = (HipoDataBank) event.getBank("REC::Event");
        HipoDataBank configBank = (HipoDataBank) event.getBank("RUN::config");
        
        helicity = eventBank.getByte("helicity", 0);
        runnum = configBank.getInt("run",0); // used for beam energy and polarization
    
        // default beam energy set to rga fall 2018
        double Eb = 10.6041; // rga fall 2018
        if (runnum > 6000) { Eb = 10.1998;} // this needs to be checked or made more precise than 6000
        // these are, at the time of writing, just for my analysis of Fall2018 inbending and Spring2019 inbending
        // I'm unsure about other beam energies at different times
        
        num_elec = recEvent.countByPid(11); // returns number of electrons
	num_piplus = recEvent.countByPid(211); 
	num_piminus = recEvent.countByPid(-211);
	num_kplus = recEvent.countByPid(321);
	num_kminus = recEvent.countByPid(-321);
        num_protons = recEvent.countByPid(2212);
        num_particles = num_elec+num_piplus+num_piminus+num_kplus+num_kminus+num_protons;
        
        // Set up Lorentz vectors
        // beam electron
        LorentzVector lv_beam = new LorentzVector();
	lv_beam.setPxPyPzM(0, 0, Math.pow(Eb*Eb-particle_mass(11)*particle_mass(11),0.5), 
                particle_mass(11));
        LorentzVector lv_target = new LorentzVector();
        // target, proton for RGA... what mass to use for RGB (deuterium target)?
	lv_target.setPxPyPzM(0, 0, 0, particle_mass(2212));
        // pull from rec banks for outgoing particles
        // electron
        String electron_index = "[11,0]"; // highest p, kinematic fitter should require FD etc
	Particle scattered_electron = recEvent.getParticle(electron_index); //
        LorentzVector lv_e = new LorentzVector();
	lv_e.setPxPyPzM(scattered_electron.px(), scattered_electron.py(), 
            scattered_electron.pz(), particle_mass(11));
        // hadrons set up below (to allow for iteration over more than two hadrons in an event)
        
        // kinematics of electron
        e_px = lv_e.px(); e_py = lv_e.py(); e_pz = lv_e.pz(); e_p = lv_e.p(); e_e = lv_e.e(); 
        e_theta = scattered_electron.theta();
        e_phi = scattered_electron.phi();
        if (e_phi < 0) { e_phi = 2*Math.PI + e_phi; }
                
        // DIS variables
        LorentzVector lv_q = new LorentzVector(lv_beam); lv_q.sub(lv_e);
	Q2 = -lv_q.mass2();
	nu = lv_beam.e()-lv_e.e();
	x  = Q2 / (2 * particle_mass(2212) * nu);
	W  = Math.pow(Math.pow(particle_mass(2212),2)+2*particle_mass(2212)*nu - Q2, 0.5);
	y = nu/lv_beam.e();
        gamma = 2*particle_mass(2212)*x/Math.pow(Q2, 0.5);
        
        // Depolarization variables
        Depolarization_A = 1/(1+gamma*gamma)*(1-y+y*y/2+y*y*gamma*gamma/4);
        Depolarization_B = 1/(1+gamma*gamma)*(1-y-y*y*gamma*gamma/4);
        Depolarization_C = (y/Math.pow(1+gamma*gamma, 0.5))*(1-y/2);
        Depolarization_V = (2-y)/(1+gamma*gamma)*Math.pow(1-y-y*y*gamma*gamma/4,0.5);
        Depolarization_W = y/(Math.pow(1+gamma*gamma, 0.5))*Math.pow(1-y-y*y*gamma*gamma/4,0.5);
    
        // set up boost to gamma*-nucleon center of mass frame
        LorentzVector gN = new LorentzVector(lv_q);
	gN.add(lv_target);
	Vector3 gNBoost = gN.boostVector();
	gNBoost.negative();
        
        // set up boost to Breit frame, this needs to be cross checked
        LorentzVector Breit = new LorentzVector(lv_q);
        LorentzVector Breit_target = new LorentzVector();
        Breit_target.setPxPyPzM(0, 0, 0, 2*x*particle_mass(2212));
        Breit.add(Breit_target);
        Vector3 BreitBoost = Breit.boostVector();
        BreitBoost.negative();
        
        // set up hadrons and dihadron
        String pIndex_string = "["+pPID+","+pIndex+"]";
        Particle hadron = recEvent.getParticle(pIndex_string);
        
        vz_e = scattered_electron.vz();
        vz_p = hadron.vz();
        
        LorentzVector lv_p = new LorentzVector();
        lv_p.setPxPyPzM(hadron.px(), hadron.py(), hadron.pz(), hadron.mass());
        
        // kinematics of hadrons
        p_px = lv_p.px(); p_py = lv_p.py(); p_pz = lv_p.pz(); p_p = lv_p.p(); p_e = hadron.e(); 
        p_theta = hadron.theta();
        p_phi = hadron.phi();
        if (p_phi < 0) { p_phi = 2*Math.PI + p_phi; }
    
        z = lv_p.e()/lv_q.e();
        
        // missing mass calculations
        LorentzVector lv_Mx = new LorentzVector(lv_q); lv_Mx.add(lv_target); lv_Mx.sub(lv_p); 
        Mx = lv_Mx.mass();
        Mx2 = lv_Mx.mass2(); // missing mass squared
        
        // boost to gamma*-nucleon center of mass frame
        LorentzVector lv_p_gN = new LorentzVector(lv_p); lv_p_gN.boost(gNBoost);
        LorentzVector lv_e_gN = new LorentzVector(lv_e); lv_e_gN.boost(gNBoost);
        LorentzVector lv_target_gN = new LorentzVector(lv_target); lv_target_gN.boost(gNBoost);
        Vector3 lv_e_gN_unit = new Vector3();
        lv_e_gN_unit.setMagThetaPhi(1, lv_e_gN.theta(), lv_e_gN.phi());
        LorentzVector lv_q_gN = new LorentzVector(lv_q); lv_q_gN.boost(gNBoost);
        Vector3 lv_q_gN_unit = new Vector3();
        lv_q_gN_unit.setMagThetaPhi(1, lv_q_gN.theta(), lv_q_gN.phi());
        Vector3 lv_y_gN = new Vector3(); // I guess this isn't really a lorenzvector so the lv name is wrong 
        lv_y_gN = lv_q_gN.vect().cross(lv_e_gN.vect()); 
        // in gamma*-nucleon frame the z axis is along gamma* and the x axis is in the 
        // e-e' plane in the direction of e. the y axis is then the cross product of these two
        
        // boost to Breit infinite momentum frame
        LorentzVector lv_p_Breit = new LorentzVector(lv_p); lv_p_Breit.boost(BreitBoost);
        LorentzVector lv_e_Breit = new LorentzVector(lv_e); lv_e_Breit.boost(BreitBoost);
        Vector3 lv_e_Breit_unit = new Vector3();
        lv_e_Breit_unit.setMagThetaPhi(1, lv_e_Breit.theta(), lv_e_Breit.phi());
        LorentzVector lv_q_Breit = new LorentzVector(lv_q); lv_q_Breit.boost(BreitBoost);
        Vector3 lv_q_Breit_unit = new Vector3();
        lv_q_Breit_unit.setMagThetaPhi(1, lv_q_Breit.theta(), lv_q_Breit.phi());
        // note that in the Breit frame +z is antialigned to the direction of q
  
        
        pT = lv_q_gN_unit.cross(lv_p_gN.vect()).mag();
    
        xF =  2*(lv_p_gN.vect() .dot(lv_q_gN.vect())) /(lv_q_gN.vect().mag()*W);
        
        zeta = lv_p_gN.e()/lv_target_gN.e(); // only really applicable when p1 is a proton
    
        p_gN_pz = lv_p_gN.vect().dot(lv_q_gN.vect())/lv_q_gN.vect().mag();
        p_Breit_pz = lv_p_Breit.vect().dot(lv_q_Breit.vect())/lv_q_Breit.vect().mag();
        
        // Breit frame rapidity
        eta = -0.5*Math.log((lv_p_Breit.e()+p_Breit_pz) /  (lv_p_Breit.e()-p_Breit_pz));
        
        // gamma*-nucleon frame rapidity
        eta_gN = 0.5*Math.log((lv_p_gN.e()+p_gN_pz) / (lv_p_gN.e()-p_gN_pz));
        
        Vector3 vecH = new Vector3();
        vecH.setMagThetaPhi(lv_p_gN.vect().mag()/z, lv_p_gN.vect().theta(), lv_p_gN.vect().phi());
        Vector3 vecR = new Vector3(vecH); // not R yet
//        vecR.setMagThetaPhi(lv_p1_gN.vect().mag()/z1, lv_p1_gN.vect().theta(), lv_p1_gN.vect().phi());
//        vecR.sub(vecH); // this is really R now that the subtraction is done
        vecR.negative();

	Vector3 R_Q = new Vector3();
							
	R_Q.setMagThetaPhi(vecR.dot(lv_q_gN_unit), lv_q_gN_unit.theta(), lv_q_gN_unit.phi());

	Vector3 vectPh = new Vector3(lv_p_gN.vect());
	Vector3 Pt_Q = new Vector3();
	Pt_Q.setMagThetaPhi(vecR.dot(lv_q_gN_unit), lv_q_gN_unit.theta(), lv_q_gN_unit.phi());
	Vector3 vectPhT = new Vector3(vectPh);
	vectPhT.sub(Pt_Q);
        
        Vector3 vT = new Vector3(lv_q_gN_unit.cross(lv_e_gN_unit)); 
        vT.unit();
	Vector3 vTH = new Vector3(lv_q_gN_unit.cross(vectPhT)); 
        vTH.unit();
        
	double cosPhiH = vT.dot(vTH);
	double sinPhiH = lv_e_gN.vect().cross(vectPhT).dot(lv_q_gN_unit);

        // scaling
	double hScale = lv_q_gN_unit.cross(lv_e_gN.vect()).mag()*lv_q_gN_unit.cross(vectPh).mag();
	sinPhiH = sinPhiH/hScale;
        
        phi = Math.acos(cosPhiH); 
        
	if (sinPhiH < 0.0) { phi = 2*Math.PI + phi; }
        
        // see trento conventions: https://arxiv.org/pdf/hep-ph/0410050.pdf        
    }
    
    
    public int get_helicity() { return -1*helicity; } // return helicity of event, 
    // -1, 0, or 1. 0 equals unassigned by EventBuilder
    
    public int get_runnum() { return runnum; }; // returns run number for polarizations and energy
    
    public int num_elec() { return num_elec; } // returns number of electrons
    
    public int num_piplus() { return num_piplus; } // returns number of piplus
    
    public int num_piminus() { return num_piminus; } // returns number of piminus
    
    public int num_kplus() { return num_kplus; }// returns number of kplus
    
    public int num_kminus() { return num_kminus; } // returns number of kminus
    
    public int num_protons() { return num_protons; } // returns number of protons
    
    public double test() { return test; } // returns test var

    public double Q2() { return Double.valueOf(Math.round(Q2*100000))/100000; } // returns Q2
    
    public double W() { return Double.valueOf(Math.round(W*100000))/100000; }// returns W
    
    public double gamma() { return Double.valueOf(Math.round(gamma*100000))/100000; } // returns gamma
    
    public double nu() { return Double.valueOf(Math.round(nu*100000))/100000; }// returns nu
    
    public double x() { return Double.valueOf(Math.round(x*100000))/100000; }// returns x
    
    public double y() { return Double.valueOf(Math.round(y*100000))/100000; }// returns y
    
    public double z() { return Double.valueOf(Math.round(z*100000))/100000; }// returns z
    
    public double Mx() { return Double.valueOf(Math.round(Mx*100000))/100000; }// returns Mx(ep1p2)
    
    public double Mx2() { return Double.valueOf(Math.round(Mx2*100000))/100000; }// returns Mx(ep1p2)
    
    public double pT() { return Double.valueOf(Math.round(pT*100000))/100000; }// returns pT
    
    public double xF() { return Double.valueOf(Math.round(xF*100000))/100000; }// returns xF
    
    public double zeta() { return Double.valueOf(Math.round(zeta*100000))/100000; }// returns zeta
    
    public double p_Breit_pz() { return Double.valueOf(Math.round(p_Breit_pz*100000))/100000; } 
    // returns p1 pz in Breit
    
    public double p_gN_pz() { return Double.valueOf(Math.round(p_gN_pz*100000))/100000; } // returns p1 pz in gN
    
    public double eta() { return Double.valueOf(Math.round(eta*100000))/100000; }// returns eta in the Breit frame
    
    public double eta_gN() { return Double.valueOf(Math.round(eta_gN*100000))/100000; }// returns eta_gN
    
    public double phi() { return Double.valueOf(Math.round(phi*100000))/100000; }
    // returns phi1 (gamma*-nucleon frame)
    
    public double Depolarization_A() { return Double.valueOf(Math.round(Depolarization_A*100000))/100000; }
    // returns Depolarization_A
    
    public double Depolarization_B() { return Double.valueOf(Math.round(Depolarization_B*100000))/100000; }
    // returns Depolarization_B
    
    public double Depolarization_C() { return Double.valueOf(Math.round(Depolarization_C*100000))/100000; }
    // returns Depolarization_C
    
    public double Depolarization_V() { return Double.valueOf(Math.round(Depolarization_V*100000))/100000; }
    // returns Depolarization_V
    
    public double Depolarization_W() { return Double.valueOf(Math.round(Depolarization_W*100000))/100000; }
    // returns Depolarization_W
    
    public double e_px() { return Double.valueOf(Math.round(e_px*100000))/100000; }// returns electron lab frame px
    
    public double e_py() { return Double.valueOf(Math.round(e_py*100000))/100000; }// returns electron lab frame py
    
    public double e_pz() { return Double.valueOf(Math.round(e_pz*100000))/100000; }// returns electron lab frame pz
    
    public double e_p() { return Double.valueOf(Math.round(e_p*100000))/100000; }// returns electron lab frame p
    
    public double e_e() { return Double.valueOf(Math.round(e_e*100000))/100000; }// returns electron lab frame energy
    
    public double e_theta() { return Double.valueOf(Math.round(e_theta*100000))/100000; } // returns electron lab 
    // frame polar angle
    
    public double e_phi() { return Double.valueOf(Math.round(e_phi*100000))/100000; } // returns electron lab 
    // frame polar angle
    
    public double p_px() { return Double.valueOf(Math.round(p_px*100000))/100000; }// returns hadron lab frame px
    
    public double p_py() { return Double.valueOf(Math.round(p_py*100000))/100000; }// returns hadron lab frame py
    
    public double p_pz() { return Double.valueOf(Math.round(p_pz*100000))/100000; }// returns hadron lab frame pz
    
    public double p_p() { return Double.valueOf(Math.round(p_p*100000))/100000; }// returns hadron lab frame p
    
    public double p_e() { return Double.valueOf(Math.round(p_e*100000))/100000; }// returns hadron lab frame energy
    
    public double p_theta() { return Double.valueOf(Math.round(p_theta*100000))/100000; } // returns p1 lab 
    // frame polar angle
    
    public double p_phi() { return Double.valueOf(Math.round(p_phi*100000))/100000; } // returns p1 lab 
    // frame polar angle
  
    public double vz_e() { return Double.valueOf(Math.round(vz_e*100000))/100000; }// returns electron z vertex
    
    public double vz_p() { return Double.valueOf(Math.round(vz_p*100000))/100000; }// returns electron z vertex
    
    
    
    private static double particle_mass (int pid) {
	if (pid==11||pid==-11) { // electron is pid=11, positron is pid=-11
            return 0.0005109989461;
	} else if (pid==211||pid==-211) { // pions
            return 0.139570;
	} else if (pid==321||pid==-321) { // kaons
            return 0.493677;
	} else if (pid==2212||pid==-2212) { // protons
            return 0.938272;
	}
            return -1;
    }
}
