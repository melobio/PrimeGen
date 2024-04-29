### NDBI

- NCBI目录： /home/x/x-multiplex-pcr/seq_search

- 代码： src

### 代码发布更新步骤

1. 更新代码

```bash
cd /home/x/x-multiplex-pcr/seq_search
// 替换 src/app.py
```

2. 构建docker

```bash
cd /home/x/x-multiplex-pcr
sudo make build-x-search
```

3. 重启 （这会重启x-search和x-pcr）

```bash
cd /home/x/x-multiplex-pcr
sudo bash run.sh
```