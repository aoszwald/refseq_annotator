import importlib.util
import sys
from pathlib import Path
from cravat import BaseAnnotator

# OC dynamically loads modules without adding their directory to sys.path.
_spec = importlib.util.spec_from_file_location(
    'refseq_transcripts_engine', Path(__file__).with_name('refseq_transcripts_engine.py'))
_engine = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_engine)


class CravatAnnotator(BaseAnnotator):
    def should_skip_chrom(self, input_data):
        return False  # Transcript cross-references do not depend on the chromosome.

    def setup(self):
        self.annotation_db = _engine.AnnotationDatabase(Path(self.data_dir) / 'mappings.sqlite')

    def annotate(self, input_data, secondary_data=None):
        return self.annotation_db.annotate(input_data)

    def cleanup(self):
        self.annotation_db.close()


if __name__ == '__main__':
    CravatAnnotator(sys.argv).run()
