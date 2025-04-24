import os
import tempfile
from drug_target_gse46602 import pipeline2
import filecmp

class TestPipeline2:

    def test_generate_probe_mappings_csv(self):
        test_resources_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "test_resources")
        probes_file = os.path.join(test_resources_dir, "probe_ids.txt")
        expected_output_file = os.path.join(test_resources_dir, "probe_mapping.csv")

        try:
            temp_dir = tempfile.mkdtemp()
            output_file = os.path.join(temp_dir, "probe_mapping.csv")

            pipeline2.generate_probe_mappings_csv(probes_file, output_file)

            assert(filecmp.cmp(expected_output_file, output_file))

        finally:
            os.remove(output_file)
            os.rmdir(temp_dir)