from pathlib import Path
from typing import Any, List, Tuple, Union

import biotite.structure as struc
import esm
import numpy as np
import torch

"""
This utils python file is almost the same as the esm_if1_utils.py that can be found in the NetTCR-struc repository.
However, small changes have been made for it to be tailored to our analysis:
- the compute_esm_if1_on_pdb function now incorporates a target_chain argument, to only run IF on a selected chain
(we only need to perform a mutational scan of the peptide)
- the forward_pass returns the full all_prob array, and not only the probability corresponding to the WT, because we 
want to perform a mutational scan on the WT structure, so we need the full picture

Some functions are kept unchanged, and some remain here even though we don't need them in our analysis 
(example: generate_esm_f1_features)
"""

def compute_esm_if1_on_pdb(
    pdb_path: Union[Path, str],
    esm_if1_model: Any,
    alphabet: Any,
    chain_names: List[str],
    target_chain: str,
) -> Tuple[List[Any], List[Any], List[str]]:
    """
    Compute ESM-IF1 encodings and probabilities on a PDB file using esm-if1 multichain.

    Args:
        pdb_path: A Pathlib object or string path to the PDB file.
        esm_if1_model: The ESM-IF1 model used for computation.
        alphabet: Alphabet for the ESM-IF1 model.
        chain_names: List of chain names in the PDB file.
        target_chain: Chain on which is performed the computation

    Returns:
        A tuple containing lists of tensor encodings, probabilities, and corresponding sequences for the PDB.
    """
    # load pdb structure into coordinates
    esm_if1_loaded_structure = esm.inverse_folding.util.load_structure(
        str(pdb_path), chain_names
    )
    (
        coords,
        native_seqs,
    ) = esm.inverse_folding.multichain_util.extract_coords_from_complex(
        esm_if1_loaded_structure
    )
    # extract esm_if1 encodings and probabilities

    esm_if1_output = {}
    seq = native_seqs[target_chain]
    esmif1_e, esmif1_p = forward_pass(
        model=esm_if1_model,
        alphabet=alphabet,
        coords=coords,
        chain=target_chain,
        seq=seq,
    )
    esm_if1_output[target_chain] = {
        "encoding": esmif1_e.cpu().detach(),
        "probability": esmif1_p.cpu().detach(),
        "sequence": seq,
    }

    return esm_if1_output


def forward_pass(
    model: Any,
    alphabet: Any,
    coords: Any,
    chain: str,
    seq: str,
):
    """
    Compute ESM-IF1 encodings and decoder probabilities for a specific chain in a PDB file.

    Args:
        model: The ESM-IF1 model used for computation.
        alphabet: Alphabet for the ESM-IF1 model.
        coords: Coordinates of the chain.
        chain: The chain identifier.
        seq: Sequence of the chain.

    Returns:
        A tuple containing the complete ESM-IF1 encoding and probability for the given chain, for all possible mutations
    """
    with torch.no_grad():
        # Preprocess coordinates
        device = next(model.parameters()).device
        all_coords = esm.inverse_folding.multichain_util._concatenate_coords(
            coords, chain
        )
        batch_converter = esm.inverse_folding.util.CoordBatchConverter(alphabet)
        batch = [(all_coords, None, seq)]
        batch_coords, confidence, _, tokens, padding_mask = batch_converter(
            batch, device=device
        )

        # Forward pass
        encoder_out = model.encoder.forward(
            batch_coords,
            padding_mask,
            confidence,
            return_all_hiddens=False,
        )
        logits, _ = model.decoder(
            tokens[:, :-1].to(device),
            encoder_out=encoder_out,
            features_only=False,
            return_all_hiddens=False,
        )

        # Prepare encoder outputs for saving
        target_chain_len = coords[chain].shape[0]
        encoder_out = encoder_out["encoder_out"][0][1:-1, 0]
        encoder_out = encoder_out[:target_chain_len].cpu().detach()

        # Prepare decoder outputs for saving
        mm = torch.nn.Softmax(dim=1)
        logits_permute = torch.squeeze(torch.permute(logits, (0, 2, 1)))

        all_probs = mm(logits_permute).cpu().detach()
        # Old way of doing it by only returning the WT probabilities
        # tgt_idxs = torch.squeeze(tokens[:, 1:]).cpu().detach()
        # res_probs = torch.stack(
        #     [all_probs[i, tgt_idxs[i].item()] for i in range(len(tgt_idxs))]
        # )
        return encoder_out, all_probs


def initialize_esm_if1_model(device) -> Tuple[Any, Any]:
    """
    Initialize the ESM-IF1 model and alphabet.

    Args:
        device: The device to use for computation.
    """
    esm_if1_model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    esm_if1_model = esm_if1_model.eval()
    esm_if1_model.to(device)
    return esm_if1_model, alphabet


def generate_esm_f1_features(
    pdb_path: Path,
    output_path: Path,
    esm_if1_model: Any,
    alphabet: Any,
    chain_names: List[str],
) -> dict:
    """Extracts structural embeddings from ESM-IF1.

    Args:
        pdb_path: Path to the PDB file.
        output_path: Path to save the extracted features.
        esm_if1_model: The ESM-IF1 model used for computation.
        alphabet: Alphabet for the ESM-IF1 model.
        chain_names: List of chain names in the PDB file.

    Returns:
        A dict with ESM-IF1 features for a structure.
    """
    # Process each row here and append the result
    esm_if1_output = compute_esm_if1_on_pdb(
        pdb_path=pdb_path,
        esm_if1_model=esm_if1_model,
        alphabet=alphabet,
        chain_names=chain_names,
    )
    torch.save(esm_if1_output, output_path)
    return esm_if1_output


def get_esm_if1_features(
    feature_path: Path,
    pdb_path: Path,
    chain_names: np.ndarray,
    esm_if1_model: Any,
    alphabet: Any,
    device: str,
) -> np.ndarray:
    """Helper function for extracting or loading ESM-IF1 features.

    Args:
        feature_path: Path to save/load features.
        pdb_path: Path to the PDB file.
        chain_id: Chain ID for each residue.
        structure: struc.AtomArray structure.
        esm_if1_model: The ESM-IF1 model used for computation.
        alphabet: Alphabet for the ESM-IF1 model.
        device: If "cuda", preprocessing is performed on the GPU.

    Returns:
        A tensor of ESM-IF1 representations.
    """
    if feature_path.exists():
        features = torch.load(feature_path, map_location=torch.device(device))
    else:
        features = generate_esm_f1_features(
            pdb_path=pdb_path,
            output_path=feature_path,
            chain_names=chain_names,
            esm_if1_model=esm_if1_model,
            alphabet=alphabet,
        )

    encodings = []
    for chain in chain_names:
        encodings.append(features[chain]["encoding"].cpu().detach())

    return torch.from_numpy(np.concatenate(encodings, axis=0)).to(device)
