#!/bin/bash
echo "remove spring pcr ..."
docker-compose -f ./mgi-pcr-adapter/docker-compose.dev.yml down -v --rmi all
echo "sucess"
