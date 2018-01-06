Building Avocado from Source
============================

You will need to have `Apache Maven <http://maven.apache.org/>`__
version 3.1.1 or later installed in order to build Avocado.

    **Note:** The default configuration is for Hadoop 2.7.3. If building
    against a different version of Hadoop, please pass
    ``-Dhadoop.version=<HADOOP_VERSION>`` to the Maven command.

.. code:: bash

    git clone https://github.com/bigdatagenomics/avocado.git
    cd avocado
    export MAVEN_OPTS="-Xmx512m -XX:MaxPermSize=128m"
    mvn clean package -DskipTests

Outputs

::

    ...
    [INFO] ------------------------------------------------------------------------
    [INFO] BUILD SUCCESS
    [INFO] ------------------------------------------------------------------------
    [INFO] Total time: 9.647s
    [INFO] Finished at: Thu May 23 15:50:42 PDT 2013
    [INFO] Final Memory: 19M/81M
    [INFO] ------------------------------------------------------------------------

You might want to take a peek at the ``scripts/jenkins-test`` script and
give it a run. We use this script to test that Avocado is
working correctly.

Running Avocado
---------------

Avocado is packaged as an
`überjar <https://maven.apache.org/plugins/maven-shade-plugin/>`__ and
includes all necessary dependencies, except for Apache Hadoop and Apache
Spark.

You might want to add the following to your ``.bashrc`` to make running
Avocado easier:

.. code:: bash

    alias avocado-submit="${AVOCADO_HOME}/bin/avocado-submit"

``$AVOCADO_HOME`` should be the path to where you have checked AVOCADO out on
your local filesystem. The alias calls a script that wraps
the ``spark-submit`` command to set up Avocado. You
will need to have the Spark binaries on your system; prebuilt binaries
can be downloaded from the `Spark
website <http://spark.apache.org/downloads.html>`__. Our `continuous
integration setup <https://amplab.cs.berkeley.edu/jenkins/job/Avocado/>`__
builds ADAM against Spark 2.0.0, Scala versions 2.10
and 2.11, and Hadoop versions 2.3.0 and 2.6.0.

Once this alias is in place, you can run Avocado by simply typing
``avocado-submit`` at the command line.

.. code:: bash

    avocado-submit

