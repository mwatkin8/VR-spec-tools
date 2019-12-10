
# VR-Spec Tool Suite

Available at: http://vcfclin.org/

The GA4GH Variant Representation Group has developed a specification ([vr-spec](https://vr-spec.readthedocs.io/en/1.0/index.html)) which includes both information models as well as an algorithm for generating globally unique identifiers. This is important because it allows for institution A to share variant information with institution B (with institution B being able to compare their VR identifiers with those provided from institution A). If they are identical then both institutions can be sure that they are talking about the exact same variant. This standard also opens the door for a variety of research and clinical tools which rely on consistent data representation.

This web suite is designed to perform the following functions:
### With a user-uploaded VCF file
    -Add VR identifiers for each variant, and return the transformed file to the user.
    -Return a the VR information models (represented in JSON) used for each entry of the VCF file.

### Without a user-uploaded VCF file
    -Return the VR information model (represented in JSON) for a user-entered HGVS expression.
    -Return the VR identifiers and the VR information model (represented in JSON) for a user-entered custom variant.
