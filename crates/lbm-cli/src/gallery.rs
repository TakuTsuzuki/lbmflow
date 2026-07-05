//! `lbm gallery` — run every built-in preset and emit a self-contained HTML
//! gallery (`index.html`, PNGs inlined as data URIs, Japanese captions).

use crate::runner;
use anyhow::{Context, Result};
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

/// Run all presets into `out_root/<preset>/` and write `out_root/index.html`.
pub fn run(out_root: &Path) -> Result<()> {
    fs::create_dir_all(out_root)
        .with_context(|| format!("出力ディレクトリを作成できません: {}", out_root.display()))?;
    let mut sections = String::new();
    for (name, desc, sc) in lbm_scenario::presets() {
        eprintln!("[gallery] {name}: 実行中…");
        let dir = out_root.join(name);
        let manifest =
            runner::run(&sc, &dir).with_context(|| format!("プリセット '{name}' の実行に失敗"))?;
        eprintln!(
            "[gallery] {name}: status={} steps={} wall={:.1}s",
            manifest.status, manifest.steps_run, manifest.wall_seconds
        );
        let mut figures = String::new();
        for f in manifest.files.iter().filter(|f| f.ends_with(".png")) {
            let bytes = fs::read(dir.join(f))
                .with_context(|| format!("PNG を読めません: {}", dir.join(f).display()))?;
            write!(
                figures,
                "<figure><img src=\"data:image/png;base64,{}\" alt=\"{}\">\
                 <figcaption>{}</figcaption></figure>\n",
                base64(&bytes),
                escape(f),
                escape(f)
            )?;
        }
        write!(
            sections,
            "<section>\n<h2>{}</h2>\n<p class=\"desc\">{}</p>\n\
             <p class=\"meta\">status={} / steps={} / {:.0} MLUPS / tau={:.3}</p>\n\
             <div class=\"figs\">\n{}</div>\n</section>\n",
            escape(name),
            escape(desc),
            escape(&manifest.status),
            manifest.steps_run,
            manifest.mlups,
            manifest.diagnostics.tau,
            figures
        )?;
    }
    let html = format!(
        r#"<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>LBMFlow プリセットギャラリー</title>
<style>
  body {{ font-family: "Hiragino Sans", "Noto Sans JP", sans-serif; margin: 2rem auto;
         max-width: 72rem; padding: 0 1rem; background: #16181d; color: #e8e8e8; }}
  h1 {{ font-size: 1.5rem; border-bottom: 2px solid #4a90d9; padding-bottom: .5rem; }}
  h2 {{ font-size: 1.15rem; color: #7cb8f2; margin-bottom: .2rem; }}
  p.desc {{ margin: .2rem 0; }}
  p.meta {{ margin: .2rem 0 .8rem; color: #9aa0aa; font-size: .85rem;
            font-family: ui-monospace, monospace; }}
  section {{ margin: 2rem 0; }}
  .figs {{ display: flex; flex-wrap: wrap; gap: 1rem; }}
  figure {{ margin: 0; }}
  figure img {{ max-width: 100%; image-rendering: pixelated; border: 1px solid #333;
                border-radius: 4px; display: block; }}
  figcaption {{ color: #9aa0aa; font-size: .8rem; font-family: ui-monospace, monospace;
                margin-top: .25rem; }}
  footer {{ margin-top: 3rem; color: #9aa0aa; font-size: .8rem; }}
</style>
</head>
<body>
<h1>LBMFlow プリセットギャラリー</h1>
<p>組み込みプリセットを <code>lbm gallery</code> で順に実行した結果のスナップショットです。
各画像は Base64 で埋め込み済み（このファイル単体で閲覧できます）。</p>
{sections}<footer>LBMFlow — 格子ボルツマン法流体シミュレータ</footer>
</body>
</html>
"#
    );
    let index = out_root.join("index.html");
    fs::write(&index, html)?;
    println!("ギャラリー生成完了: {}", index.display());
    Ok(())
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

/// Minimal RFC 4648 base64 (standard alphabet, padded); avoids a dependency.
fn base64(data: &[u8]) -> String {
    const TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
    for chunk in data.chunks(3) {
        let n = (u32::from(chunk[0]) << 16)
            | (u32::from(*chunk.get(1).unwrap_or(&0)) << 8)
            | u32::from(*chunk.get(2).unwrap_or(&0));
        out.push(TABLE[(n >> 18) as usize & 63] as char);
        out.push(TABLE[(n >> 12) as usize & 63] as char);
        out.push(if chunk.len() > 1 {
            TABLE[(n >> 6) as usize & 63] as char
        } else {
            '='
        });
        out.push(if chunk.len() > 2 {
            TABLE[n as usize & 63] as char
        } else {
            '='
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base64_matches_rfc4648_vectors() {
        // RFC 4648 test vectors
        for (input, expect) in [
            (&b""[..], ""),
            (b"f", "Zg=="),
            (b"fo", "Zm8="),
            (b"foo", "Zm9v"),
            (b"foob", "Zm9vYg=="),
            (b"fooba", "Zm9vYmE="),
            (b"foobar", "Zm9vYmFy"),
        ] {
            assert_eq!(base64(input), expect);
        }
    }

    #[test]
    fn escape_covers_html_metacharacters() {
        assert_eq!(escape(r#"a<b>&"c""#), "a&lt;b&gt;&amp;&quot;c&quot;");
    }
}
