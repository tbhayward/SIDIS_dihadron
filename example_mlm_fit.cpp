// Timothy B. Hayward MLM fit
// calculate asymmetry for pTpT for back-to-back proton-pi+ analysis
#include <iostream>
#include <iomanip>
#include <fstream>

int current_bin = 0;

const int num_bins = 6;
float counts[6];
float means[6];

// 1-D Quantile Bins

float pTpT_bins[7] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6}; // bin edges (so only six A_LU points here)

TTree *tree = new TTree("T","B2B yields");

int b2b_counts, status, runnum, evnum, helicity;
float e_p, e_theta, e_phi, vz_e;
float pi_p, pi_theta, pi_phi, vz_pi;
float P_p, P_theta, P_phi, vz_P;
float Q2, W, x, y, z_pi, z_P;
float Mx, Mx1, Mx2;
float zeta, Mh;
float pT_pi, pT_P, pTpT;
float xF_pi, xF_P, eta_pi, eta_P, Delta_eta;
float phi_pi, phi_P, Delta_phi;

float pdf(int runnum, float Delta_phi, double *par) {

	float pol; 
	if (runnum == 11 ) { pol = 0.86; } // MC
	else if (runnum >= 5032 && runnum < 5333) { pol = 0.8592; } 
	else if (runnum >= 5333 && runnum <= 5666) { pol = 0.8922; }
	else if (runnum >= 6616 && runnum <= 6783) { pol = 0.8453; }
	// polarization is added in event-by-event 
	// (so we don't divide by an average polarization at the end)

	// B2B analysis fit
	return pol*(par[0]*TMath::Sin(Delta_phi)+par[1]*TMath::Sin(2*Delta_phi));
}

// negative log-likelihood function
void negative_log_likelihood(Int_t &npar, Double_t *gin, Double_t &nll, Double_t *par,
	Int_t iflag) {

	int num_pairs = tree->GetEntries(); // total number of hadron pairs (basically events, but
	// we could have more than one pair of hadrons per trigger electron)

	// separate sums for both helicity states
	double sum_N = 0;
	double sum_P = 0;

	for (int i=0; i<num_pairs; i++) {
	// for (int i=0; i<25000; i++) { // if I 
		if (i % 100000 == 0 ) { 
			cout << "On bin: "<<current_bin << " with " <<  " ";
			cout << sum_N << " " << sum_P << ". ";}
			// this is kind of pointless, I don't really know what the numbers mean during 
			// intermediate steps but it's fun to watch them change and converge on the final value
		if (i==0) {
			cout << endl << endl << endl << endl;
			cout << current_bin << " ";
		}
		tree->GetEntry(i);

		float Delta_Y = Delta_eta; // difference in rapidities but some places I used eta and some
		// places I used Y so this is necessary

		// B2B cuts (w/ Harut Avakian)
		// my event cuts (I process a more general sample of the data)
		if (helicity == 0 ) { continue; }
		if (status == 1e0) { continue; } // check for QA

		if (e_theta < 10*3.14159/180 || e_theta > 30*3.14159/180) { continue; }
		if (pi_theta > 30*3.14159/180) { continue; }
		if (P_theta > 30*3.14159/180) { continue; }
		
		if (y > 0.75) { continue; }
		if (pi_p < 1.20 || pi_p > 4.0) { continue; }
		if (P_p < 0.5) { continue; }

		// vertex cuts
		if (TMath::Abs(vz_e + 3) > 5) { continue; }
		if (TMath::Abs(vz_e - vz_pi) > 7) { continue; }
		if (TMath::Abs(vz_e - vz_P) > 7) { continue; }
		if (TMath::Abs(vz_pi - vz_P) > 5) { continue; }

		if (z_pi < 0.2 || z_pi > 0.7) { continue; }

		// ensure back-to-back events, these cuts are hotly contested (particularly by Miguel)
		// may drop
		if (xF_P > 0.0) { continue; }
		if (xF_pi < 0.0) { continue; }
		if (Delta_eta < 0.00) { continue; }

		// cut out exclusive pi- and rho- events
		if (Mx < 0.85) { continue; }
		if (Mx2 < 1.6) { continue; }

		// determine bin
		if (pTpT < pTpT_bins[current_bin] || pTpT > pTpT_bins[current_bin+1]) { continue; }

		// this is the actual building of the likelihood
		// make sure you get the helicity signs correct, the skims have the helicity backwards (I think)
		// but when I write my text files I reverse the sign
		// can always check the sign of the SIDIS pi+ asymmetry (should be positive)
		if (helicity == -1) {  // helicity definition is flipped
			sum_N += TMath::Log(1-pdf(runnum, Delta_phi, par));
		} else if (helicity == 1) {
			sum_P += TMath::Log(1+pdf(runnum, Delta_phi, par));
		}
	}

	// negative log likelihood, normalization minus sum of both states
	nll =  num_pairs - sum_P - sum_N;
}

void mlm_fits() {

	float bins[num_bins+1] = {}; bins[0] = 0; 
	for (int i = 0; i < num_bins; ++i) { 
		bins[i+1] = pTpT_bins[i+1]; 
	}

	float results[12][num_bins][3] = {{{}}};

	// FALL 2018 // 
	// TString filename = "/volatile/clas12/thayward/scratch/rga_fall2018_inbending_proton_piplus.txt";
	// TString filename_original = "/work/clas12/thayward/SIDIS/yields/proton_piplus_analysis/proton_piplus/rga_fall2018_inbending_proton_piplus.txt";

	// SPRING 2019 // 
	// TString filename = "/volatile/clas12/thayward/scratch/rga_spring2019_inbending_proton_piplus.txt";
	// TString filename_original = "/work/clas12/thayward/SIDIS/yields/proton_piplus_analysis/proton_piplus/rga_spring2019_inbending_proton_piplus.txt";

	// COMBINED DATA // 

	TString filename = "/volatile/clas12/thayward/scratch/rga_combined_inbending_proton_piplus_g.txt";
	TString filename_original = "/work/clas12/thayward/SIDIS/yields/proton_piplus_analysis/proton_piplus/rga_combined_inbending_proton_piplus.txt";


	TString command = "cp "+filename_original+" "+filename;
	system(command);

	FILE *fp = fopen(filename,"r");
	TFile *hfile = hfile = TFile::Open(filename,"RECREATE");
	// TTree *tree = new TTree("T","dihadron yields");
	tree->Branch("status",&status,"status/I"); tree->SetBranchAddress("status", &status);
	tree->Branch("runnum",&runnum,"runnum/I"); tree->SetBranchAddress("runnum", &runnum);
	tree->Branch("evnum",&evnum,"evnum/I"); tree->SetBranchAddress("evnum", &evnum);
  	tree->Branch("helicity",&helicity,"helicity/I"); tree->SetBranchAddress("helicity", &helicity);
  	tree->Branch("e_p",&e_p,"e_p/F"); tree->SetBranchAddress("e_p", &e_p);
  	tree->Branch("e_theta",&e_theta,"e_theta/F"); tree->SetBranchAddress("e_theta", &e_theta);
  	tree->Branch("e_phi",&e_phi,"e_phi/F"); tree->SetBranchAddress("e_phi", &e_phi);
  	tree->Branch("vz_e",&vz_e,"vz_e/F"); tree->SetBranchAddress("vz_e", &vz_e);
  	tree->Branch("pi_p",&pi_p,"pi_p/F"); tree->SetBranchAddress("pi_p", &pi_p);
  	tree->Branch("pi_theta",&pi_theta,"pi_theta/F"); tree->SetBranchAddress("pi_theta", &pi_theta);
  	tree->Branch("pi_phi",&pi_phi,"pi_phi/F"); tree->SetBranchAddress("pi_phi", &pi_phi);
  	tree->Branch("vz_pi",&vz_pi,"vz_pi/F"); tree->SetBranchAddress("vz_pi", &vz_pi);
  	tree->Branch("P_p",&P_p,"P_p/F"); tree->SetBranchAddress("P_p", &P_p);
  	tree->Branch("P_theta",&P_theta,"P_theta/F"); tree->SetBranchAddress("P_theta", &P_theta);
  	tree->Branch("P_phi",&P_phi,"P_phi/F"); tree->SetBranchAddress("P_phi", &P_phi);
  	tree->Branch("vz_P",&vz_P,"vz_P/F"); tree->SetBranchAddress("vz_P", &vz_P);
  	tree->Branch("Q2",&Q2,"Q2/F"); tree->SetBranchAddress("Q2", &Q2);
  	tree->Branch("W",&W,"W/F"); tree->SetBranchAddress("W", &W);
  	tree->Branch("x",&x,"x/F"); tree->SetBranchAddress("x", &x);
  	tree->Branch("y",&y,"y/F"); tree->SetBranchAddress("y", &y);
  	tree->Branch("z_pi",&z_pi,"z_pi/F"); tree->SetBranchAddress("z_pi", &z_pi);
  	tree->Branch("z_P",&z_P,"z_P/F"); tree->SetBranchAddress("z_P", &z_P);
  	tree->Branch("Mx",&Mx,"Mx/F"); tree->SetBranchAddress("Mx", &Mx);
  	tree->Branch("Mx1",&Mx1,"Mx1/F"); tree->SetBranchAddress("Mx1", &Mx1);
  	tree->Branch("Mx2",&Mx2,"Mx2/F"); tree->SetBranchAddress("Mx2", &Mx2);
  	tree->Branch("zeta",&zeta,"zeta/F"); tree->SetBranchAddress("zeta", &zeta);
  	tree->Branch("Mh",&Mh,"Mh/F"); tree->SetBranchAddress("Mh", &Mh);
  	tree->Branch("pT_pi",&pT_pi,"pT_pi/F"); tree->SetBranchAddress("pT_pi", &pT_pi);
  	tree->Branch("pT_P",&pT_P,"pT_P/F"); tree->SetBranchAddress("pT_P", &pT_P);
  	tree->Branch("pTpT",&pTpT,"pTpT/F"); tree->SetBranchAddress("pTpT", &pTpT);
  	tree->Branch("xF_pi",&xF_pi,"xF_pi/F"); tree->SetBranchAddress("xF_pi", &xF_pi);
  	tree->Branch("xF_P",&xF_P,"xF_P/F"); tree->SetBranchAddress("xF_P", &xF_P);
  	tree->Branch("eta_pi",&eta_pi,"eta_pi/F"); tree->SetBranchAddress("eta_pi", &eta_pi);
  	tree->Branch("eta_P",&eta_P,"eta_P/F"); tree->SetBranchAddress("eta_P", &eta_P);
  	tree->Branch("Delta_eta",&Delta_eta,"Delta_eta/F"); tree->SetBranchAddress("Delta_eta", &Delta_eta);
  	tree->Branch("phi_pi",&phi_pi,"phi_pi/F"); tree->SetBranchAddress("phi_pi", &phi_pi);
  	tree->Branch("phi_P",&phi_P,"phi_P/F"); tree->SetBranchAddress("phi_P", &phi_P);
  	tree->Branch("Delta_phi",&Delta_phi,"Delta_phi/F"); tree->SetBranchAddress("Delta_phi", &Delta_phi);
  	

  	
  	char line[450];
  	int i = 0;
	while (fgets(line,450,fp)) {
		sscanf(&line[0],"%d %d %d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e ",
		&status, &runnum, &evnum, &helicity,
		&e_p, &e_theta, &e_phi, &vz_e,
		&pi_p, &pi_theta, &pi_phi, &vz_pi,
		&P_p, &P_theta, &P_phi, &vz_P,
		&Q2, &W, &x, &y, &z_pi, &z_P,
		&Mx, &Mx1, &Mx2, &zeta, &Mh,
		&pT_pi, &pT_P, &pTpT,
		&xF_pi, &xF_P, &eta_pi, &eta_P, &Delta_eta,
		&phi_pi, &phi_P, &Delta_phi);
		tree->Fill();
	}
	tree->Write();


	int nPar = 2; // number of fit parameters, here I fit A\sin\Delta\phi + B\sin(2\Delta\phi)

	// mlm fit
	for (int i = 0; i < num_bins; i++) {
		bool check = false;
		for (int j=0; j<tree->GetEntries(); j++) {

			tree->GetEntry(j);

			// B2B cuts (w/ Harut Avakian)
			if (helicity == 0 ) { continue; }
			if (status == 1e0) { continue; } // check for QA

			if (e_theta < 10*3.14159/180 || e_theta > 30*3.14159/180) { continue; }
			if (pi_theta > 30*3.14159/180) { continue; }
			if (P_theta > 30*3.14159/180) { continue; }
			
			if (y > 0.75) { continue; }
			if (pi_p < 1.20 || pi_p > 4.0) { continue; }
			if (P_p < 0.5) { continue; }

			if (TMath::Abs(vz_e + 3) > 5) { continue; }
			if (TMath::Abs(vz_e - vz_pi) > 7) { continue; }
			if (TMath::Abs(vz_e - vz_P) > 7) { continue; }
			if (TMath::Abs(vz_pi - vz_P) > 5) { continue; }

			if (z_pi < 0.2 || z_pi > 0.7) { continue; }

			if (xF_P > 0.0) { continue; }
			if (xF_pi < 0.0) { continue; }
			if (Delta_eta < 0.00) { continue; }

			if (Mx < 0.85) { continue; }
			if (Mx2 < 1.6) { continue; }

			// find correct bin
			if (pTpT < pTpT_bins[current_bin] || pTpT > pTpT_bins[current_bin+1]) { continue; }

			// I'm calculating the mean value of the bins to place the asymmetry points at
			counts[current_bin]++;
			if (check) {
				means[current_bin]=means[current_bin]+pTpT;
			} else {
				means[current_bin]=pTpT;
				check = true;
			}
		}

		// prepare minuit
		double arglist[10]; arglist[0] = 1; 
   		int ierflg = 0;
		TMinuit m(nPar);
		m.SetFCN(negative_log_likelihood);
		m.SetPrintLevel(0); // -1 quiet, 0 normal, 1 verbose
		// 1 for chi2 fit, 0.5 for negative log-likelihood fit
		// see section 1.4.1 in MINUIT manual, e.g., http://hep.fi.infn.it/minuit.pdf
		m.SetErrorDef(0.5); 
		// VERY IMPORTANT, USE 1 for chi2 fit and 0.5 for MLM fit
		// !!!!!!!!!!!!!!!!!!!!!!!

		// parameter no., name, start value, step size, range min., range max.
		// range min = range max = 0 -> no limits
		m.DefineParameter(0, "par0", -0.020, 0.01, -0.2, 0.1);
		m.DefineParameter(1, "par1", 0.000, 0.01, -0.3, 0.3);
		// make sure your results 
		arglist[0] = 2000;      // max calls

		// // now ready for minimization step
		m.Migrad();
		m.mnexcm("MINImize",arglist,1,ierflg);
		m.mnexcm("MINOS",arglist,1,ierflg);
		m.mnexcm("HESSE",arglist,0,ierflg);
		m.Command("SHOW COV"); // show covariance matrix

		double A0, A0err;
		double A1, A1err;

		m.GetParameter(0, A0, A0err);
		m.GetParameter(1, A1, A1err);

		results[0][i][1] = A0; results[0][i][2] = A0err;
		results[1][i][1] = A1; results[1][i][2] = A1err;
		current_bin++;

		cout << endl << endl << endl << endl << endl;

		// print out the results
		cout << "{";
		for (int j = 0; j < nPar; ++j) {
			cout << "{";
			for (int k = 0; k < num_bins; ++k) {
				if (k < num_bins-1) {
					cout << "{" << means[k]/counts[k] << ", " << results[j][k][1] << ", " << 
						results[j][k][2] << "}, ";
				} else {
					cout << "{" << means[k]/counts[k] << ", " << results[j][k][1] << ", " << 
						results[j][k][2] << "}";
				}
			}
			if (j < nPar-1) { cout << "},"; } else {cout << "}";}
		}
		cout << "};";


		cout << endl << endl << endl << endl;
		for (int j = 0; j < nPar; ++j) {
			for (int k = 0; k < num_bins; ++k) {
				cout << counts[k] << " " << means[k]/counts[k] << ", " << results[j][k][1] << ", " << 
					results[j][k][2] << endl;
			}
			cout << endl;
		}

		cout << endl << endl << endl << endl << endl << endl << endl << endl;

	}

}