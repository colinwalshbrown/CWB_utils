#!/usr/bin/env python
"""
module to handle annotation features (may or may not be tied to a sequence / alignment)
"""

class AnnotationTrack:
    """
    set of annotations for a single contiguous piece of sequence
    """
    def __init__(self,length=0,features=None,parentSeq=None):
        self.length = length
        self.features = features
        self.parentSeq = parentSeq
        self.tags = {}

    def addFeature(self,feature):
        if not self.features:
            self.features = [feature]
        else:
            self.features.append(feature)

    def addTag(self,tag_key_value_pair):
        self.tags[tag_key_value_pair[0]] = tag_key_value_pair[1]

    def setParentSeq(self,sequence):
        if (len(sequence) == self.length):
            self.parentSeq = sequence
        else:
            raise AnnotationError("Sequence length does not match annotation length")

    def getAllFeatures(self):
        if self.features:
            return self.features
        else:
            return []

    def getFeaturesByType(self,matchType):
        type_matches = []
        if self.features:
            type_matches = [x for x in self.features if x.type == matchType]
        return type_matches
    
    def getFeaturesInRegion(self,matchStart,matchEnd):
        region_matches = [x for x in self.features if ((x.start >= matchStart) and (x.end <=matchEnd))]
        return region_matches
    
    def newTrackFromType(self,matchType):
        matches = self.getFeaturesByType(matchType)
        return AnnotationTrack(length=self.seqlength,features=matches,parentSeq=self.parentSeq)
    
    def newTrackFromRegion(self,matchStart,matchEnd):
        matches = self.getFeaturesInRegion(matchStart,matchEnd)
        return AnnotationTrack(length=(matchEnd - matchStart),features=matches,parentSeq=self.parSeq)
        
    def newTrackFromScore(self,scoreCutoff):
        matches = [x for x in self.features if x.score >= scoreCutoff] 

    def getMaxFeature(self,tag):
        return max(self.features,key=(lambda x: x.tags[tag]))

class _Feature:
    """
    Base class for features
    """
    def __init__(self,start=0,end=0,type="",tags=None):
        self.length = (end - start)
        self.start = start
        self.end = end
        self.type = type
        self.tags = tags
    
    def addTag(self,tag):
        if not self.tags:
            self.tags = {tag[0]:tag[1]}
        else:
            self.tags[tag[0]] = tag[1]
            
    def __str__(self):
        rep = "%s|start:%s|end:%s|tags:" % (self.type,self.start,self.end)
        tag_str = ",".join(["(" + str(x[0]) + ":" + str(x[1]) + ")" for x in self.tags.items()])
        return rep + tag_str

class GenericFeature(_Feature):
    def __init__(self):
        _Feature._init__(self,type="GenericFeature")
        
class BindingSiteFeature(_Feature):
    def __init__(self,start=0,end=0,type="BindingSiteFeature",tags={}):
        _Feature.__init__(self,start,end,type,tags)

class GeneFeature(_Feature):
    def __init__(self):
        _Feature.__init__(self,type="GeneFeature")

class ExonFeature(_Feature):
    def __init__(self):
        _Feature.__init__(self,type="ExonFeature")

class IntronFeature(_Feature):
    def __init__(self):
        _Feature.__init__(self,type="IntronFeature")

class AnnotationError(Exception):
    def __init__(self, value):
        self.val = value
    def __str__(self):
        return repr(self.val)

        
    
    
    
    
    
    
    
