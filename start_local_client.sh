#!/bin/bash
echo "start local client"
source scripts/export-env.sh .env;
source scripts/export-env.sh .env.local;
yarn start:client
echo "success"

