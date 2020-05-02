# Copyright (c) 2013-2020, SIB - Swiss Institute of Bioinformatics and
#                          Biozentrum - University of Basel
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#   http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


"""
Test example codes from documentation.
Each example code is given as a python script in scripts-folder and may use
data from data-folder.
Scripts are executed with ost within BUILDFOLDER/tests/doc.
Output of scripts shall be checked and cleaned up by caller.
make target: test_doctests.py_run
Don't forget to add new data-files and scripts to CMakeLists.txt!
"""
import unittest
import os
import subprocess
import shutil
from ost import io, settings, table
import qmean

class DocTests(unittest.TestCase):

    def runScript(self, script_path, arguments=[]):
        """Run script with python binary and given arguments.
        Returns tuple (return-code, stdout, stderr).
        This script assumes that all required Python libraries (ost, qmean)
        can be found in default locations or in the properly set 
        PYTHON_PATH environment variable. If you perform a default unit
        test run (i.e. make check), this is taken care off by cmake 
        (see qmean_unittest macro).
        Additionally, cmake sets the PYTHON_BINARY environment variable
        with which we execute scripts. If PYTHON_BINARY cannot be found,
        the default python from the system is used as a fallback (In the hope
        that ost and qmean are compiled with that one).
        AND YES, CALLING IT WITH RAW PYTHON MEANS NO COMPOUND LIBRARY ETC.
        """

        python_bin = os.getenv("QMEAN_PYTHON_BINARY")
        if python_bin is None:
            raise RuntimeError("Must set QMEAN_PYTHON_BINARY env variable "
                               "to run unit tests. CMake should take care "
                               "of that when executing make check.")

        # launch it
        cmd = [python_bin, script_path] + arguments
        job = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        sout, serr = job.communicate()
        return job.returncode, sout.decode(), serr.decode()


    def testTorsionPotentialExample(self):
        # run it 
        return_code, sout, serr = self.runScript('torsion_potential_example.py')

        # check whether script ran through fine and generated the requested files
        self.assertEqual(return_code, 0)
        self.assertTrue(os.path.exists("torsion_stat.dat"))
        self.assertTrue(os.path.exists("torsion_pot.dat"))

        # load the generated files and check the options
        loaded_pot = qmean.TorsionPotential.Load("torsion_pot.dat")
        loaded_stat = qmean.TorsionStatistic.Load("torsion_stat.dat")

        pot_group_identifier = loaded_pot.GetOptions().group_identifier
        stat_group_identifier = loaded_stat.GetOptions().group_identifier
        expected_group_identifier = ["all-GLY-all", "all-PRO-all", 
                                     "all-all-PRO", "all-VAL,ILE-all", 
                                     "all-all-all"]

        for a,b,c in zip(expected_group_identifier, pot_group_identifier, 
                         stat_group_identifier):
            self.assertTrue(a == b and a == c)

        pot_bins = loaded_pot.GetOptions().number_of_bins
        stat_bins = loaded_stat.GetOptions().number_of_bins
        expected_bins = [1,1,18,18,1,1]

        for a,b,c in zip(expected_bins, pot_bins, stat_bins):
            self.assertTrue(a == b and a == c)

        # check some characteristics specific to crambin
        self.assertEqual(loaded_stat.GetTotalCount(), 46-4) 
        self.assertEqual(loaded_stat.GetTotalCount("all-GLY-all"), 4)
        self.assertEqual(loaded_stat.GetTotalCount("all-PRO-all"), 5)
        self.assertEqual(loaded_stat.GetTotalCount("all-all-PRO"), 5)
        self.assertEqual(loaded_stat.GetTotalCount("all-VAL,ILE-all"), 6)
        self.assertEqual(loaded_stat.GetTotalCount("all-all-all"), 22)

        # cleanup
        os.remove("torsion_stat.dat")
        os.remove("torsion_pot.dat")

    def testSSAgreementExample(self):
        # run it 
        return_code, sout, serr = self.runScript('ss_agreement_example.py')
        self.assertEqual(return_code, 0)

        # lets parse the values we get
        sout = sout.splitlines()
        values = list()
        for line in sout:
            if "Score given a DSSP state of" in line:
                values.append(float(line.split()[-1]))

        expected_values = [0.107144355774, -1.07871556282, 0.136189222336, 
                           -1.01866531372, 0.224040657282, -1.0713198185, 
                           0.302175998688, -1.25567805767, 0.412762403488, 
                           -1.38902068138, 0.541625022888, -1.72601318359, 
                           0.666330933571, -2.01082897186, 0.797769427299, 
                           -2.49906778336, 0.92312681675, -3.09363245964, 
                           1.09601712227, -5.09709453583]

        self.assertEqual(len(values), len(expected_values))

        for a,b in zip(values, expected_values):
            self.assertAlmostEqual(a,b,4)

    def testSmootherExample(self):
        return_code, sout, serr = self.runScript('smoother_example.py')
        self.assertEqual(return_code, 0)

    def testPSIPREDHandlerExample(self):
        return_code, sout, serr = self.runScript('psipred_handler_example.py')
        self.assertEqual(return_code, 0)
        sout = sout.splitlines()

        # braindead reading of the standard out and comparing to what it once was
        # on my machine. Values can change if SSAgreement scorer is retrained, feel
        # free to adapt the unit test in that case.
        full_chain_pred = sout[1]
        full_chain_conf = sout[2]
        subset_pred = sout[4]
        subset_conf= sout[5]
        full_chain_scores = [float(item) for item in sout[7].strip('[]').split(',')]
        subset_scores = [float(item) for item in sout[8].strip('[]').split(',')]

        exp_full_chain_pred = "CCCCCCHHHHHHHHHEECCCCCHHHHHHHCCCEECCCCCCCCCCCC"
        exp_full_chain_conf = "9657578887665410368988878977449676559978188898"
        exp_subset_pred = "CHHHHHHHHHEECC"
        exp_subset_conf = "78887665410368"
        exp_full_chain_scores = [0.7959088683128357, -0.7178850173950195, 
        -0.5056777000427246, 0.7294700741767883, 0.5967527627944946, 
        0.8243091702461243, 0.9231268167495728, 0.9231268167495728, 
        0.9231268167495728, 0.7977694272994995, 0.6663309335708618, 
        0.6663309335708618, 0.5416250228881836, 0.4127624034881592, 
        0.13618922233581543, -1.536504864692688, -2.0208404064178467, 
        0.5437932014465332, 0.6971776485443115, 0.9875132441520691, 
        0.7852041721343994, 0.7852041721343994, 0.9231268167495728, 
        0.7977694272994995, 0.9231268167495728, 1.0960171222686768, 
        0.7977694272994995, 0.7977694272994995, 0.4127624034881592, 
        -1.1485421657562256, 0.8720449209213257, 0.656849205493927, 
        1.2513508796691895, 1.1609821319580078, 0.5699062943458557, 
        0.5967527627944946, 0.8720449209213257, 0.8720449209213257, 
        0.7294700741767883, 0.7852041721343994, 0.26415103673934937, 
        -0.3378578722476959, -0.3378578722476959, -0.3378578722476959, 
        0.7959088683128357, 0.7852041721343994]
        exp_subset_scores = [0.8243091702461243, 0.9231268167495728, 
        0.9231268167495728, 0.9231268167495728, 0.7977694272994995, 
        0.6663309335708618, 0.6663309335708618, 0.5416250228881836, 
        0.4127624034881592, 0.13618922233581543, -1.536504864692688, 
        -2.0208404064178467, 0.5437932014465332, 0.6971776485443115]

        self.assertEqual(full_chain_pred, exp_full_chain_pred)
        self.assertEqual(full_chain_conf, exp_full_chain_conf)
        self.assertEqual(subset_pred, exp_subset_pred)
        self.assertEqual(subset_conf, exp_subset_conf)
        self.assertEqual(len(full_chain_scores), len(exp_full_chain_scores))
        self.assertEqual(len(subset_scores), len(exp_subset_scores))

        for a,b in zip(full_chain_scores, exp_full_chain_scores):
            self.assertAlmostEqual(a,b,4)

        for a,b in zip(subset_scores, exp_subset_scores):
            self.assertAlmostEqual(a,b,4)

    def testPackingPotentialExample(self):
        return_code, sout, serr = self.runScript('packing_potential_example.py')
        self.assertEqual(return_code, 0)

        # check for generated files
        self.assertTrue(os.path.exists("packing_statistic.dat"))
        self.assertTrue(os.path.exists("packing_potential.dat"))
        self.assertTrue(os.path.exists("packing_energy_plot.png"))

        stat = qmean.PackingStatistic.Load("packing_statistic.dat")
        pot = qmean.PackingPotential.Load("packing_potential.dat")

        self.assertAlmostEqual(stat.GetOptions().cutoff, pot.GetOptions().cutoff, 4)
        self.assertEqual(stat.GetOptions().max_counts, pot.GetOptions().max_counts)
        self.assertAlmostEqual(pot.GetOptions().sigma, 0.02, 4)

        ca_pseudo_energies = list()
        nz_pseudo_energies = list()

        for count in range(33):
           e = pot.GetEnergy(qmean.ChemType.C_LYS_A, count)
           ca_pseudo_energies.append(e)
           e = pot.GetEnergy(qmean.ChemType.N_LYS_Z, count)
           nz_pseudo_energies.append(e)

        # for low counts we expect z_pseudo_energies to be lower, for
        # high counts we expect it to be the other way arount.
        self.assertTrue(ca_pseudo_energies[0] > nz_pseudo_energies[0])
        self.assertTrue(ca_pseudo_energies[-1] < nz_pseudo_energies[-1])

        os.remove("packing_statistic.dat")
        os.remove("packing_potential.dat")
        os.remove("packing_energy_plot.png")

    def testMembraneExample(self):
        return_code, sout, serr = self.runScript('membrane_example.py')
        self.assertEqual(return_code, 0)
        self.assertTrue(os.path.exists("gpcr.pdb"))
        self.assertTrue(os.path.exists("gpcr_transmembrane_part.pdb"))      
        self.assertTrue(os.path.exists("energy_gap_plot.png"))  
        os.remove("gpcr.pdb")
        os.remove("gpcr_transmembrane_part.pdb")
        os.remove("energy_gap_plot.png") 

    def testLocalScorerExample(self):
        return_code, sout, serr = self.runScript('local_scorer_example.py')
        self.assertEqual(return_code, 0)
        self.assertTrue(os.path.exists("local_scorer.dat"))
        os.remove("local_scorer.dat")

    def testInteractionPotentialExample(self):
        return_code, sout, serr = self.runScript('interaction_potential_example.py')
        self.assertEqual(return_code, 0)
        self.assertTrue(os.path.exists("interaction_statistic.dat"))
        self.assertTrue(os.path.exists("interaction_potential.dat"))
        self.assertTrue(os.path.exists("interaction_energy_plot.png"))

        pot = qmean.InteractionPotential.Load("interaction_potential.dat")
        stat = qmean.InteractionStatistic.Load("interaction_statistic.dat")

        self.assertAlmostEqual(stat.GetOptions().lower_cutoff, pot.GetOptions().lower_cutoff, 4)
        self.assertAlmostEqual(stat.GetOptions().upper_cutoff, pot.GetOptions().upper_cutoff, 4)
        self.assertEqual(stat.GetOptions().number_of_bins, pot.GetOptions().number_of_bins)
        self.assertEqual(stat.GetOptions().sequence_sep, pot.GetOptions().sequence_sep)
        self.assertAlmostEqual(pot.GetOptions().sigma, 0.02, 4)

        # brain dead parsing of standard out... just check whether the crambin overall energy 
        # matches with the expected value
        crambin_e = float(sout.splitlines()[1].split()[-1])
        exp_crambin_e = -0.0329360701144
        self.assertAlmostEqual(crambin_e, exp_crambin_e, 4)

        os.remove("interaction_statistic.dat")
        os.remove("interaction_potential.dat")
        os.remove("interaction_energy_plot.png")

    def testDiscoExample(self):
        return_code, sout, serr = self.runScript('disco_example.py')
        self.assertEqual(return_code, 0)
        self.assertTrue(os.path.exists("assigned_disco.pdb"))
        self.assertTrue(os.path.exists("distance_constraint.png"))
        self.assertTrue(os.path.exists("dc.dat"))

        exp_scores = [0.821259319782,0.783250153065,0.772677242756,0.759754657745,
        0.784342288971,0.806196153164,0.792514741421,0.779994547367,0.837361574173,
        0.826949357986,0.780048549175,0.836500585079,0.829025328159,0.771560072899,
        0.745820641518,0.428177326918,0.421403765678,0.300698012114,0.375659704208,
        0.249514341354,0.792763710022,0.837131500244,0.823297321796,0.786858141422,
        0.83382666111,0.878767430782,0.876249074936,0.867045164108,0.849434673786,
        0.792085945606,0.775034964085,0.755913853645,0.770276725292,0.680759072304,
        0.826633691788,0.682561337948,0.721345543861,0.830392122269,0.833950400352,
        0.856562495232,0.864427745342,0.87565600872,0.857636272907,0.845507740974,
        0.826569795609,0.817613482475,0.813940405846,0.767517626286,0.793221354485,
        0.793132781982,0.833293199539,0.816828370094,0.82264560461,0.842022359371,
        0.76438587904,0.692659318447,0.539588809013,0.268566578627,0.31766885519,
        0.342374056578,0.422271996737,0.387684702873,0.44497436285,0.502450227737,
        0.581216990948,0.525600552559,0.599538445473,0.687871396542,0.750949621201,
        0.845186948776,0.85482609272,0.717804610729,0.792412221432,0.864990830421,
        0.801212012768,0.818390727043,0.848885774612,0.874608099461,0.837935864925,
        0.863742411137,0.883148491383,0.889209568501,0.863461196423,0.833316743374,
        0.85692769289,0.839185655117,0.83184915781,0.845653712749,0.880329668522,
        0.849484324455,0.865452587605,0.834379673004,0.845095276833,0.835442602634,
        0.863957166672,0.835496068001,0.863817155361,0.880394816399,0.869758486748,
        0.84107708931,0.833675205708,0.887885212898,0.874950945377,0.870262384415,
        0.87125223875,0.865008175373,0.871275484562,0.881080389023,0.878927409649,
        0.836000323296,0.844562888145,0.869691371918,0.872847139835,0.865396618843,
        0.867342233658,0.876452982426,0.865706682205,0.858490586281,0.870484113693,
        0.865277290344,0.867626845837,0.859246730804,0.87319457531,0.870645642281,
        0.857725322247,0.851319789886,0.83846527338,0.811100900173,0.82908308506,
        0.828980624676,0.816096127033,0.76682394743,0.724892377853,0.806805789471,
        0.837520897388,0.840397536755,0.818115890026,0.858565032482,0.848163306713,
        0.837713003159,0.848052501678,0.850918412209,0.862653255463,0.87482625246,
        0.862853288651,0.839833557606,0.844416320324,0.835634112358,0.844650566578,
        0.798519432545,0.80936384201,0.839978575706,0.832882404327,0.750095963478,
        0.777846932411,0.794464051723,0.795528411865,0.820773005486,0.811093568802,
        0.82362729311,0.816468238831,0.847517192364,0.863434195518,0.860754609108,
        0.828978955746,0.794631004333,0.767912626266,0.828493893147,0.85519272089,
        0.8408203125,0.805510699749,0.818524599075,0.843751370907,0.822210609913,
        0.80524623394,0.78876721859,0.749441146851,0.625212907791,0.513365626335,
        0.479973077774,0.740260720253,0.863315999508,0.859677374363,0.838380515575,
        0.837450683117,0.837235152721,0.835918903351,0.812751412392,0.846859633923,
        0.859999358654,0.874414086342,0.815671443939,0.845019817352,0.868504583836,
        0.872323393822,0.865004479885,0.886147916317,0.873253107071,0.883003890514,
        0.865707755089,0.783871471882,0.656572759151,0.481571108103,0.590996146202,
        0.433227926493,0.344746947289,0.347757399082,0.437933683395,0.440915167332,
        0.517445445061,0.623940289021,0.595308482647,0.667969286442,0.671737551689,
        0.650144279003,0.794988751411,0.820125639439,0.815876245499,0.809833645821,
        0.816134214401,0.815550327301,0.837843954563,0.852966427803,0.860906362534,
        0.854808807373,0.815634250641,0.817847371101,0.843975186348,0.859192311764,
        0.848325192928,0.840126156807,0.76463419199,0.79028236866,0.847461223602,
        0.830110609531,0.884908437729,0.88754093647,0.889272987843,0.879056930542,
        0.86722022295,0.827437698841,0.844816744328,0.855136394501,0.849663555622,
        0.831016302109,0.824648976326,0.857223331928,0.889174699783,0.884983003139,
        0.896271586418,0.887469291687,0.869247972965,0.857521474361,0.853989601135,
        0.84389680624,0.859813988209,0.836564898491,0.858998835087,0.885324597359,
        0.895750403404,0.898643255234,0.884318292141,0.891741514206,0.878203749657,
        0.870936870575,0.873240172863,0.876082956791,0.868219673634,0.824413537979,
        0.815245747566,0.814845860004,0.854783892632,0.829991936684,0.813286304474,
        0.84471321106,0.849817812443,0.805973947048,0.847839415073,0.834975719452]

        model = io.LoadPDB("assigned_disco.pdb")
        self.assertEqual(len(exp_scores), len(model.residues))
        for exp_score, r in zip(exp_scores, model.residues):
            # storing value as befactor limits accuracy
            self.assertAlmostEqual(exp_score, r.atoms[0].b_factor, 2)

        # load the distance constraints that we saved to disk, redo scoring and
        # check again
        dc = qmean.DisCoContainer.Load("dc.dat")
        scores = dc.GetScores(model.CreateFullView())
        self.assertEqual(len(scores), len(exp_scores))
        for exp_score, score in zip(exp_scores, scores):
            # lossy compression when saving DisCoContainer limits accuracy
            self.assertAlmostEqual(exp_score, score, 2)

        os.remove("assigned_disco.pdb")
        os.remove("distance_constraint.png")
        os.remove("dc.dat")

    def testContainerExample(self):
        return_code, sout, serr = self.runScript('container_example.py')
        self.assertEqual(return_code, 0)

        self.assertTrue(os.path.exists("potential_container.dat"))
        os.remove("potential_container.dat")

        # in the example script we generate a container with potentials, estimate
        # some energies, dump the container, load it again and repeat the energy
        # calculation with the loaded container (which should be exactly the same)
        # let's check the expected equality.
        sout = sout.splitlines()
        int_e_before = float(sout[0].split()[-1])
        pack_e_before = float(sout[1].split()[-1])
        int_e_after = float(sout[4].split()[-1])
        pack_e_after = float(sout[5].split()[-1])
        self.assertEqual(int_e_before, int_e_after)
        self.assertEqual(pack_e_before, pack_e_after)

    def testChemtypeExample(self):
        return_code, sout, serr = self.runScript('chemtype_example.py')
        self.assertEqual(return_code, 0)

    def testCBetaPotentialExample(self):
        return_code, sout, serr = self.runScript('cbeta_potential_example.py')
        self.assertEqual(return_code, 0)
        # we don't test too much here as its similar to the interaction potential
        self.assertTrue(os.path.exists("cbeta_energy_plot.png"))
        os.remove("cbeta_energy_plot.png")

    def testACCPROHandlerExample(self):
        return_code, sout, serr = self.runScript('accpro_handler_example.py')
        self.assertEqual(return_code, 0)
        sout = sout.splitlines()

        # braindead reading of the standard out and comparing to what it once was
        # on my machine. Values can change if SSAgreement scorer is retrained, feel
        # free to adapt the unit test in that case.
        full_chain_pred = sout[1]
        subset_pred = sout[3]
        full_chain_scores = [float(item) for item in sout[5].strip('[]').split(',')]
        subset_scores = [float(item) for item in sout[6].strip('[]').split(',')]

        exp_full_chain_pred = "bbbbbbeeeeeebbbbbebebbebebebbebbbbbbbbbeeeeebb"
        exp_subset_pred = "beeeeeebbbbbeb"
        exp_full_chain_scores = [0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 
        1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
        exp_subset_scores = [1.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 
        1.0, 0.0, 1.0, 0.0]

        self.assertEqual(full_chain_pred, exp_full_chain_pred)
        self.assertEqual(subset_pred, exp_subset_pred)
        self.assertEqual(len(full_chain_scores), len(exp_full_chain_scores))
        self.assertEqual(len(subset_scores), len(exp_subset_scores))

        for a,b in zip(full_chain_scores, exp_full_chain_scores):
            self.assertAlmostEqual(a,b,4)

        for a,b in zip(subset_scores, exp_subset_scores):
            self.assertAlmostEqual(a,b,4)

    def testAssessModelQualityExample(self):
        return_code, sout, serr = self.runScript('assess_model_quality_example.py')
        self.assertEqual(return_code, 0)

        self.assertTrue(os.path.exists("local_profile.png"))
        self.assertTrue(os.path.exists("qmean4_ref_plot.png"))
        self.assertTrue(os.path.exists("qmean6_ref_plot.png"))
        self.assertTrue(os.path.exists("qmean4_sliders.png"))
        self.assertTrue(os.path.exists("qmean6_sliders.png"))
        os.remove("local_profile.png")
        os.remove("qmean4_ref_plot.png")
        os.remove("qmean6_ref_plot.png")
        os.remove("qmean4_sliders.png")
        os.remove("qmean6_sliders.png")

        qmean4_score = float(sout.splitlines()[0].split()[-1])
        qmean6_score = float(sout.splitlines()[1].split()[-1])
        qmean4_z_score = float(sout.splitlines()[2].split()[-1])
        qmean6_z_score = float(sout.splitlines()[3].split()[-1])
        avg_local_score = float(sout.splitlines()[5].split()[-1])
        avg_local_score_error = float(sout.splitlines()[6].split()[-1])
        exp_qmean4_score = 0.765673453311
        exp_qmean6_score = 0.755831800470
        exp_qmean4_z_score = -0.1425450733423469
        exp_qmean6_z_score = -0.25876529067689186
        exp_avg_local_score = 0.866393307804
        self.assertAlmostEqual(qmean4_score, exp_qmean4_score, 4)
        self.assertAlmostEqual(qmean6_score, exp_qmean6_score, 4)
        self.assertAlmostEqual(qmean4_z_score, exp_qmean4_z_score, 4)
        self.assertAlmostEqual(qmean6_z_score, exp_qmean6_z_score, 4)
        self.assertAlmostEqual(avg_local_score, exp_avg_local_score, 4)
        # check for NaN
        self.assertFalse(avg_local_score_error == avg_local_score_error) 

        local_score_line = sout.splitlines()[4]
        local_scores = eval(local_score_line[local_score_line.index('{'):])
        exp_local_scores = {'A': {1: 0.7617946984290086, 2: 0.8103660238896163, 
        3: 0.9159270528774698, 4: 0.9471315272930895, 5: 0.9219615800683376, 
        6: 0.907399880118726, 7: 0.8542986822209679, 8: 0.9291329364368498, 
        9: 0.9903208136662552, 10: 0.8635250330267045, 11: 0.9280658100174551, 
        12: 0.9293427248017916, 13: 0.9288973329440364, 14: 0.8848821357205228, 
        15: 0.8798075825212647, 16: 0.890420110918607, 17: 0.7737988119699544, 
        18: 0.7821949181451636, 19: 0.813352159068194, 20: 0.8966737492003395, 
        21: 0.8794009556650662, 22: 0.9005752723031046, 23: 0.863339873844473, 
        24: 0.9145065919166018, 25: 0.8298285290543451, 26: 0.9096768373406201, 
        27: 0.9640837138716855, 28: 0.8679302119718261, 29: 0.8564178798145665, 
        30: 0.9487899692621875, 31: 0.9811979389469101, 32: 0.9736665885010267, 
        33: 0.9116909302976822, 34: 0.8722555291551453, 35: 0.8220538336900243, 
        36: 0.7995659798522442, 37: 0.7931476732106149, 38: 0.7952036724601584, 
        39: 0.7573134943937226, 40: 0.8364165103492325, 41: 0.8136952379960155, 
        42: 0.7806667545942363, 43: 0.8086675483501075, 44: 0.8212366973904793, 
        45: 0.7547270179166898, 46: 0.7887433534989572}}
        self.assertTrue('A' in local_scores)
        for k,v in exp_local_scores['A'].items():
            self.assertTrue(k in local_scores['A'])
            self.assertAlmostEqual(local_scores['A'][k], v, 2)


    def testAssessMembraneModelQualityExample(self):
        return_code, sout, serr = self.runScript('assess_membrane_model_quality.py')
        self.assertEqual(return_code, 0)

        self.assertTrue(os.path.isdir('original_hhblits_alignment'))
        self.assertTrue(os.path.isdir('shift_in_front_helix_four'))
        self.assertTrue(os.path.isdir('shift_into_middle'))
        self.assertTrue(os.path.isdir('shift_towards_cter'))

        return_code, sout, serr = self.runScript('reproduce_fig4_from_publication.py')
        self.assertEqual(return_code, 0)
        self.assertTrue(os.path.exists("alignment_comparison.png"))

        def _CheckResult(name):

            tab_path = os.path.join("example_data", name + "_local_scores.txt")
            exp_tab = table.Table.Load(tab_path)
            tab_path = os.path.join(name, "local_scores.txt")
            tab = table.Table.Load(tab_path)
            exp_qmean_col = exp_tab["QMEAN"] 
            qmean_col = tab["QMEAN"] 
            self.assertEqual(len(exp_qmean_col), len(qmean_col))

            # we accept a really low number of discrepancies
            # the reason is that the membrane finding can have slight
            # differences leading to residues that are suddenly 
            # classified differently (membrane/interface/soluble). This results
            # in different statical potentials being applied and therefore 
            # different numbers...
            n_fails = 0
            for a,b in zip(exp_qmean_col, qmean_col):
                if abs(a-b) > 0.01:
                    n_fails += 1
            self.assertTrue(n_fails<=3)

        _CheckResult('original_hhblits_alignment')
        _CheckResult('shift_in_front_helix_four')
        _CheckResult('shift_into_middle')
        _CheckResult('shift_towards_cter')
        
        shutil.rmtree('original_hhblits_alignment')
        shutil.rmtree('shift_in_front_helix_four')
        shutil.rmtree('shift_into_middle')
        shutil.rmtree('shift_towards_cter')
        os.remove('alignment_comparison.png')

    def testRegressorTraining(self):
        return_code, sout, serr = self.runScript('regressor_training.py')
        if return_code != 0:
            if "ModuleNotFoundError" in serr:
                print("Could not import keras/pandas, skip regressor training test")
                return
        self.assertEqual(return_code, 0)
        testing_rmse = float(sout.splitlines()[-1].split()[-1])
        # there is some randomness in training but we should for sure have an 
        # rmse below 10 if everything went well...
        self.assertTrue(testing_rmse < 10.0)
         
if __name__ == "__main__":
    from ost import testutils
    testutils.RunTests()
