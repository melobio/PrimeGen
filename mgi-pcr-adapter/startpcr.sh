#!/bin/bash
echo "start spring pcr ...."
mvn -f ./mgi-pcr-adapter/pom.xml install -Dmaven.skip.test=true
mvn -f ./mgi-pcr-adapter/cloud-pacs-primer/pom.xml install -Dmaven.skip.test=true
docker-compose -f ./mgi-pcr-adapter/docker-compose.dev.yml up -d
echo "sucess"