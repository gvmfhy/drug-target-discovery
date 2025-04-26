import os
import tempfile
from drug_target_gse46602 import pipeline2
import filecmp

class TestPipeline2:

    test_resources_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_resources")

    def test_generate_probe_mappings_csv(self):
        # Given an input file containing a list of probe ids
        probes_file = os.path.join(self.test_resources_dir, "tiny_probe_ids.txt")
        # And given an output file to write to
        temp_dir = tempfile.mkdtemp()
        output_file = os.path.join(temp_dir, "probe_mapping.csv")

        try:
            # When I attempt to generate a mappings file from those probes
            pipeline2.generate_probe_mappings_csv(probes_file, output_file)

            # Then the designated output file should match what we have recorded from a previous run.
            expected_output_file = os.path.join(self.test_resources_dir, "tiny_probe_mapping.csv")
            assert(filecmp.cmp(expected_output_file, output_file))

        finally:
            # Cleanup temp files
            os.remove(output_file)
            os.rmdir(temp_dir)

    def test_extract_probe_mapping_from_file(self):
        # Given I have a list of probe ids and a corresponding csv file mapping those ids
        probe_ids = ["1007_s_at", "1053_at", "117_at", "121_at", "1255_g_at", "1294_at", "1316_at"]

        # When I attempt to extract those mappings from the csv file
        probes_csv = os.path.join(self.test_resources_dir, "tiny_probe_mapping.csv")
        mapping = pipeline2.extract_probe_mapping_from_file(probe_ids, probes_csv)

        # Then I get a mapping that looks like this:
        expected_mapping = {"1007_s_at": "DDR1", 
                            "1053_at":"RFC2",
                            "117_at":"HSPA6",
                            "121_at":"PAX8",
                            "1255_g_at":"GUCA1A",
                            "1294_at":"UBA7",
                            "1316_at":"THRA"}
        assert expected_mapping == mapping