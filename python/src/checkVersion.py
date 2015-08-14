from openturns import *

def isAfter(neededMajor, neededMinor, neededRevision):
    revision = '0'
    listVersion = PlatformInfo.GetVersion().split('.')
    major = listVersion[0]
    minor = listVersion[1]
    if len(listVersion) == 3:
        revision = listVersion[2]
        
    # Parse the major
    majorInt = -1
    try:
        majorInt = int(major)
    except:
        print "Warning! Strange major number", major
        return False
    if majorInt < neededMajor:
        return False
    # Parse the minor
    minorInt = -1
    try:
        minorInt = int(minor)
    except:
        print "Warning! Strange minor number", minor
        return False
    if minorInt < neededMinor:
        return False
    # Parse the revision
    revisionInt = -1
    try:
        revisionInt = int(revision)
    except:
        # Should be a devel revision, try to extract the associated official revision
        try:
            # Promote the revision to the next one
            revisionParts = revision.split('-')
            revisionInt = int(revisionParts[0]) + 1
        except:
            print "Warning! Strange revision number", revision
            return False
    if revisionInt <= neededRevision:
        return False
    return True
