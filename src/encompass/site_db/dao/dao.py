# dao.py

from sqlalchemy.orm import Session
from encompass.site_db.models.models import ChainStructure, PdbChainInformation, PdbStructureWhole

# ============================================================
# CREATE
# ============================================================

def save_pdbwhole(db: Session, pdb_obj: PdbStructureWhole):
    """Insert a PDB object (with children) into DB."""
    db.add(pdb_obj)
    db.commit()
    db.refresh(pdb_obj)
    return pdb_obj


# ============================================================
# READ
# ============================================================

def get_wholepdbByPDBCode(db: Session, pdb: str):
    return db.query(PdbStructureWhole).filter(PdbStructureWhole.pdb_code == pdb).first()

# ============================================================
# READ
# ============================================================

def get_wholepdb(db: Session, psw_id: int):
    return db.query(PdbStructureWhole).filter(PdbStructureWhole.psw_id == psw_id).first()


def get_all_wholepdbs(db: Session):
    return db.query(PdbStructureWhole).all()


# ============================================================
# UPDATE
# ============================================================

def update_wholepdb(db: Session, pdb_obj: PdbStructureWhole):
    """Update using a PDB object already modified."""
    db.commit()
    db.refresh(pdb_obj)
    return pdb_obj


# ============================================================
# DELETE
# ============================================================

def delete_pdb(db: Session, pdb_obj: PdbStructureWhole):
    db.delete(pdb_obj)
    db.commit()
    return True

# ============================================================
# CREATE PDB Chain
# ============================================================

def save_pdbChain(db: Session, chain_obj: ChainStructure):
    """Insert a PDB object (with children) into DB."""
    db.add(chain_obj)
    db.commit()
    db.refresh(chain_obj)
    return chain_obj
