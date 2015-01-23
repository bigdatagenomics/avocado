# avocado Internals

This section discusses the hooks provided in avocado for extending avocado. For users
interested in the details of avocado specific algorithms, please see [Algorithm Details].

## Extending avocado

avocado provides hooks for creating new read pre-processing, variant discovery, genotyping,
and variant post-processing stages. To create a new version of a stage, the user must
define the following:

* *Preprocessing:* An object that extends the
`org.bdgenomics.avocado.preprocessing.PreprocessingStage` trait
* *Variant discovery:* An object that extends the
`org.bdgenomics.avocado.discovery.ExplorerCompanion` trait, and a class that extends
the `org.bdgenomics.avocado.discovery.Explorer` trait.
* *Genotyping:* An object that extends the
`org.bdgenomics.avocado.genotyping.GenotyperCompanion` trait, and a class that extends
the `org.bdgenomics.avocado.genotyping.Genotyper` trait.
* *Postprocessing:* An object that extends the
`org.bdgenomics.avocado.postprocessing.PostprocessingStage` trait, and a class that
extends the `org.bdgenomics.avocado.postprocessing.GenotypeFilter` trait.

The object will be passed a `HierarchicalConfiguration` that corresponds to the
instance of the stage being executed. From here, the stage should extract configuration
parameters.

Currently, we require users to manually modify the list of stages that exist. Long term,
we will identify stages by inspecting the classpath. This will make it easier for users
to create proprietary stages.

## Statistics

The `AvocadoConfigAndStats` class is meant to provide access to global values that may
be generally useful. All of the values computed in this class are declared as lazy; this
means that they are only computed when used. This choice was made because avocado is designed
as a configurable pipeline. If the values were not lazy, this would force them to be computed
on class initialization, which could be expensive if the pipeline that had been instantiated
did _not_ make use of values in the statistics class.

A single instance of the `AvocadoConfigAndStats` class is created as soon as reads are loaded
in. Once this instance is created, it is passed to every single stage initialization.
