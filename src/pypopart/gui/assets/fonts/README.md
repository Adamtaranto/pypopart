# Bundled fonts

Self-hosted so the GUI keeps its typography offline — it runs as a local
server, and a CDN font would silently fall back to a system face on a
machine with no network.

| File                         | Family                    | Used for                         | Source                                                                      |
| ---------------------------- | ------------------------- | -------------------------------- | --------------------------------------------------------------------------- |
| `ArchivoBlack-Regular.woff2` | Archivo Black             | App title                        | [google/fonts](https://github.com/google/fonts/tree/main/ofl/archivoblack)  |
| `SpaceGrotesk.woff2`         | Space Grotesk (variable)  | UI text, network labels          | [google/fonts](https://github.com/google/fonts/tree/main/ofl/spacegrotesk)  |
| `JetBrainsMono.woff2`        | JetBrains Mono (variable) | Alignment viewer, mutation ticks | [google/fonts](https://github.com/google/fonts/tree/main/ofl/jetbrainsmono) |

All three are licensed under the SIL Open Font License 1.1; see `OFL.txt`.
Converted from the upstream TTFs to woff2 with `fontTools`.
