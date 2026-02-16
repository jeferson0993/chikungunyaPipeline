# Base Image: Python 3.9 Slim (Debian-based)
FROM python:3.9-slim

# Metadados do Container
LABEL maintainer="Jeferson F Silva <jeferson0993@gmail.com>"
LABEL version="1.0"
LABEL description="Pipeline Genômico para Chikungunya com MAFFT e Biopython"

# 1. Instalar dependências de sistema (MAFFT é crucial aqui)
# O flag -y confirma automaticamente e rm -rf limpa o cache para reduzir tamanho da imagem
RUN apt-get update && apt-get install -y \
    mafft \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# 2. Configurar diretório de trabalho
WORKDIR /app

# 3. Copiar e instalar dependências Python
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 4. Copiar o código fonte
COPY src/ ./src/

# 5. Criar diretório para volumes (persistência de dados)
RUN mkdir -p /app/data

# 6. Definir variáveis de ambiente para evitar arquivos .pyc e buffer de log
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# 7. Comando padrão de execução
CMD ["python", "src/main.py"]
