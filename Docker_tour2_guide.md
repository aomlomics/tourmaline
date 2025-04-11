## Docker

To load the image on a computer with Docker installed:

```
docker load < tourmaline2-dev-amd64.tar.gz 
```

To see the image:

```
docker images
```

Run the local image in a container:

```
docker run -v $HOME:/data -it tourmaline2-dev-amd64
```

## 