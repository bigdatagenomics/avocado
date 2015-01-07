#!/bin/sh

mvn -Dresume=false -P distribution release:clean release:prepare release:perform
