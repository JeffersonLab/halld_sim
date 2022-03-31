<pre>
   Author:        Curtis A. Meyer

   Creation Date: March 2022

   NOTE:    This code is written in python3, and will only work properly if the HDDM python
            wrappers are built for python3 and are added to the PYTHONPATH


   Purpose: This code will generate photoproduction events of the form:

                           gamma p -> p X  -> p a b c

            where the intermediate X will decay to three final state mesons.
            If any of a,b and c are pi0 or eta mesons, those will be decayed
            to two photons.

            The photproduction reaction chooses a photon energy uniformly distributed
            between 8 GeV and 9 GeV, and selects the mass of X uniformly between
            1.02(mass(a) + mass(b) + mass(c)) and 3.0GeV. The photoproduction reaction
            is then produced with a fixed t-slope, where the default is 2, but this can
            be changed with the -t or --tslope options.

            The decay to the three photon final state is thrown according to three body
            phase space, such that for a fixed mass(X), the Dalitz plots will be uniformly
            populated.

            The program currently has a number of possibiliies for a,b and c which are
            selected using the -r or --reaction options. The current ones are

            -r pi0pi0pi0 generates pi0 pi0 pi0, where all three pi0s decay to 2 photons.
            -r pi+pi-pi0 generates pi+ pi- pi0, where the pi0 decays to 2 photons.
            -r pi+pi-eta generates pi+ pi- eta, where the eta decays to 2 photons.
            -r pi0pi0eta generates pi0 pi0 eta, where the eta and pi0s all decay to 2 photons.
            -r etaetapi0 generates eta eta pi0, where the eta's and pi0 all decay to 2 photons.
            -r etaetaeta generates eta eta eta, where the eta's all decay to 2 photons.
            -r k+k-pi0 generates K+ K- pi0, where the pi0 decays to 2 photons.

            Other command-line options control the run number, number of events and output
            file name.
            -R or --run RunNo sets the run number, the default is 99999.
            -e or --events NoEvt sets the numbre of events to generate, default is 1,000,000.
            -f or --file FileName sets the output file name.

   Caveats: This code is written in python3, so your PYTHONPATH system variable must contain
            the directory  "${HALLS_RECON_HOME}/${BMS_OSNAME}/python3"

            _________________________________________________________________________________
</pre>

