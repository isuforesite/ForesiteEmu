# ForesiteEmu

Latin Hypercube Sampling for computer experiment emulation as part of the 
FORESITE project.

`ForesiteEmu` is an R package that generates samples to emulate a crop computer 
model as a combination of Latin Hypercube Samples from three family of inputs: 
soil variables, weather variables, management decisions. 
This software is part of the FORESITE efforts.

### Install

```
devtools::install_github("luisdamiano/ForesiteEmu", build_vignettes = TRUE)
```

### Prerequisites
  * R (>= 3.5.0)
  
### Examples

See the vignette:

```
vignette("example", package = "ForesiteEmu")
```
  
## Contributing

Please, reach out to maintainer to coordinate efforts.

## Authors

* **Luis Damiano** - *Creator, author, maintainer* - [luisdamiano](https://github.com/luisdamiano)

## License
GPL (>=3)

## Acknowledgments

* We acknowledge support from Iowa State University’s Presidential
Interdisciplinary Research Initiative.