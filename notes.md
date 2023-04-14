
Originally wrote when these config files all lived in `~/src/metabolike`. Need to test
docker up process doesn't depend on other stuff in that directory.

To start Docker container (making if it doesn't exist yet):
```
sudo docker compose up -d
```

To create database, after Docker command above:
```
conda activate metabolike
metabolike setup -c metabolike-setup.yaml --no-create-db
```
Note that despite the `--no-create-db` flag, which seemed necessary, it does seem to be
making a database...


To clean up Docker data (to remake from scratch):
```
sudo docker compose down -v
```


