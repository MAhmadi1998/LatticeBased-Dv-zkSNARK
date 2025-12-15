# 📜 Lattice-Based Designated Verifier zkSNARK

This repository contains the source code for the paper:

**“Lattice-Based Designated Verifier zkSNARKs From Standard Assumptions”**
https://eprint.iacr.org/2025/2099

> This work is based on the implementation from  
> [module_gaussian_lattice](https://github.com/lucasprabel/module_gaussian_lattice), which provides foundational support for lattice trapdoors on modules and their applications.

---

## 🔧 Getting Started

### 🔨 Compilation

To compile the main protocol file, simply run:

```bash
make mainprotocol.c
```

### ▶️ Running the Protocol

After compiling, execute the protocol with:

```bash
./mainprotocol
```

---

## ⚙️ Configuration

You can **customize the length of the witness** by changing the following parameter in the `common.h` file:

```c
#define PARAM_NumTrapdoors <desired_length>
```

Replace `<desired_length>` with your desired witness length.

---

## 📁 Repository Structure

```
.
├── common.h              # Shared parameters (e.g., PARAM_NumTrapdoors)
├── mainprotocol.c        # Main protocol source file
├── Makefile              # Build configuration
└── ...
```

---

## 🧠 Acknowledgments

- Based on [lucasprabel/module_gaussian_lattice](https://github.com/lucasprabel/module_gaussian_lattice)
- Thanks to the cryptographic community for continued research in post-quantum secure zkSNARKs.
