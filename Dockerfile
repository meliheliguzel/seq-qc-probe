# ==========================================
# Builder Stage (Build and Dependencies)
# ==========================================
FROM python:3.13-slim-bookworm AS builder

COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

WORKDIR /app

COPY pyproject.toml uv.lock ./

ENV UV_COMPILE_BYTECODE=1
ENV UV_LINK_MODE=copy
RUN uv sync --frozen --no-install-project --no-dev

COPY src/ ./src
COPY README.md ./
RUN uv pip install --no-deps .

# ==========================================
# Runner Stage (Lightweight Runtime)
# ==========================================
FROM python:3.13-slim-bookworm AS runner

LABEL org.opencontainers.image.title="seq-qc-probe" \
      org.opencontainers.image.description="High-throughput streaming pre-flight FASTQ quality control gate" \
      org.opencontainers.image.url="https://github.com/meliheliguzel/seq-qc-probe" \
      org.opencontainers.image.source="https://github.com/meliheliguzel/seq-qc-probe" \
      org.opencontainers.image.licenses="MIT"

RUN groupadd -r bio && useradd -r -g bio -u 1000 bio

WORKDIR /app

COPY --from=builder --chown=bio:bio /app/.venv /app/.venv

ENV PATH="/app/.venv/bin:$PATH"
ENV PYTHONUNBUFFERED=1

USER bio
WORKDIR /data

ENTRYPOINT ["seq-qc-probe"]
CMD ["--help"]
