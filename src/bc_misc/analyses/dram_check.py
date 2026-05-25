# %%
"""
Check on the DRAM

Format as in
https://git.psu.edu/BlackCAT/shared/interfaces/-/blob/main/dram-pixel-format.md

"""
import matplotlib as mpl
mpl.use('module://ipympl.backend_nbagg')

from tools import *
import gzip

dramdir = top_dir.joinpath("payload_checkout/REDOWNLOADED_DRAMS/dramApr29nd")


def readsplit(dir, letters: str):
    data = []
    for c in letters:
        f = next(dir.glob(f"*{c}.S"))
        data.append(np.frombuffer(gzip.open(f, "rb").read(), dtype=np.uint32))
    return np.concatenate(data)


def hexes(
    d,
    item_format="{:08x}",
    linelength=0,
    sep=" ",
    linesep="\n",
    incipit="",
    terminus="\n",
):
    strs = [item_format.format(v) for v in d]
    if linelength > 0:
        chunks = [strs[i : i + linelength] for i in range(0, len(strs), linelength)]
    else:
        chunks = [strs]
    result = incipit + linesep.join([sep.join(chunk) for chunk in chunks]) + terminus
    return result

def field(structs:list[dict], fieldname:str|list[str]):
    if isinstance(fieldname, str):
        return [struct[fieldname] for struct in structs]
    return [field(structs, name) for name in fieldname]
    
    

# %%
dA = readsplit(dramdir, "A")
print(hexes(dA[0:10], linelength=4))

# %%
# header tags are in pix.rs


def word_to_range(w):
    return ((w & 0x3FF), (w >> 10) & 0x3FF)


def frame_header(data, istart=0, start_address=0):
    words = data[istart:]
    tags = (words[0:] >> 24).astype(np.uint8)
    if tags[0] != 0xFE:
        raise RuntimeError("Frame not starting with 0xfe tag")
    # What the next tag will be
    # 0=normal, 1=back pointer 2=forward pointer
    pointer_state = 0
    for idx, (tag, word) in enumerate(zip(tags, words)):
        match pointer_state:
            case 0:
                pass
            case 1:
                back_pointer = word
                pointer_state = 0
                continue
            case 2:
                forward_pointer = word
                pointer_state = 0
                continue
        match tag:
            case 0xFE:
                match (word >> 17) & 0b0111:
                    case 0b000:
                        n_hist_words = 0
                    case 0b111:
                        n_hist_words = 1024
                    case _:
                        raise RuntimeError("Wrong histogram size")

                is_dense = (word & 0x1_0000) != 0
                write_epoch = word & 0xFFFF
                pointer_state = 1
            case 0xFC:
                n_pixel_words = word & 0x00FF_FFFF
                pointer_state = 2
            case 0xC0:
                frame_idx = (word >> 4) & 0xF_FFF
                frame_status = word & 0xF
            case 0xC1:
                timestamp0 = word & 0xFF_FFFF
            case 0xC2:
                timestamp1 = word & 0xFF_FFFF
            case 0xC3:
                timestamp2 = word & 0x00_FFFF
            case 0xC7:
                correlation_tag = word & 0xFF
            case 0xC8:
                hist_rows = word_to_range(word)
            case 0xC9:
                hist_cols = word_to_range(word)
            case 0xCA:
                pixel_rows = word_to_range(word)
            case 0xCB:
                pixel_cols = word_to_range(word)
            case 0xC5:
                break
            case _:
                raise RuntimeError(f"Unexpected Tag at {idx} {tag:02x}")
    else:
        raise RuntimeError(f"Did not find end of header after {idx} words")
    timestamp = (timestamp2 << 48) | (timestamp1 << 24) | timestamp0
    result = dict(
        timestamp=timestamp,
        back_pointer=back_pointer,
        forward_pointer=forward_pointer,
        n_hist_words=n_hist_words,
        n_pixel_words=n_pixel_words,
        is_dense=is_dense,
        write_epoch=write_epoch,
        correlation_tag=correlation_tag,
        hist_rows=hist_rows,
        hist_cols=hist_cols,
        pixel_rows=pixel_rows,
        pixel_cols=pixel_cols,
        frame_idx=frame_idx,
        frame_status=frame_status,
    )
    return result


def headers(data, start_address=0):
    tags = (data >> 24).astype(np.uint8)
    w = np.ravel(np.argwhere(tags == 0xFE))
    headers = [frame_header(data, istart) for istart in w]
    return headers

heads = headers(data=dA)
times = np.array(field(heads,'timestamp'))*1e-8

for h in heads:
    print(h)


plt.plot(times, field(heads, 'frame_idx'), '.')

# %%
