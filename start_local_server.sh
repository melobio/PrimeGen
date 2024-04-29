#!/bin/bash
echo "start local server"
source scripts/export-env.sh .env;
source scripts/export-env.sh .env.local;
yarn start:server
echo "success"
