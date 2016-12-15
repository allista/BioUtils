from Bio.SeqFeature import SeqFeature

class AnnotatorBase(object):
    annotation_type = 'NONE'
    feature_type = 'misc_feature'
    tag_qualifier = 'annotation_tag'
    type_qualifier = 'annotation_type'
    len_quaifier = 'annotation_length'

    def annotate_location(self, name, group, location):
        feature = SeqFeature(location, type=self.feature_type)
        feature.qualifiers['ugene_group'] = group
        feature.qualifiers['ugene_name'] = name
        feature.qualifiers[self.tag_qualifier] = name
        feature.qualifiers[self.type_qualifier] = self.annotation_type
        feature.qualifiers[self.len_quaifier] = len(location)
        return feature

    @classmethod
    def check_type(cls, feature, atype=None):
        try: t = feature.qualifiers[cls.type_qualifier][0]
        except (KeyError, IndexError): return False
        return t == (atype or cls.annotation_type)

    @classmethod
    def get_tag(cls, feature):
        try: tag = feature.qualifiers[cls.tag_qualifier][0]
        except (KeyError, IndexError): tag = None
        return tag

class HSPAnnotator(AnnotatorBase):
    score_qualifier = 'annotation_score'

    def hsp_score(self, hsp): return 0

    def hsp2feature(self, name, group, location, hsp):
        feature = self.annotate_location(name, group, location)
        feature.qualifiers[self.score_qualifier] = self.hsp_score(hsp)
        return feature

    @classmethod
    def get_score(cls, feature):
        try: score = float(feature.qualifiers[cls.score_qualifier][0])
        except (KeyError, IndexError, ValueError): score = 0
        return score
