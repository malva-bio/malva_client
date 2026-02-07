"""Tests for malva_client.config."""

import json
import os
import pytest

from malva_client.config import Config


class TestConfigLoad:
    def test_load_from_file(self, tmp_path):
        config_file = tmp_path / "config.json"
        config_file.write_text(json.dumps({
            "server_url": "https://malva.example.com",
            "api_token": "tok-123",
            "verify_ssl": False,
        }))

        config = Config.load(str(config_file))
        assert config.server_url == "https://malva.example.com"
        assert config.api_token == "tok-123"
        assert config.verify_ssl is False

    def test_load_from_env_vars(self, tmp_path, monkeypatch):
        # Use a non-existent config path so only env vars are used
        config_file = tmp_path / "nope.json"
        monkeypatch.setenv("MALVA_API_URL", "https://env.malva.bio")
        monkeypatch.setenv("MALVA_API_TOKEN", "env-token")

        config = Config.load(str(config_file))
        assert config.server_url == "https://env.malva.bio"
        assert config.api_token == "env-token"

    def test_env_overrides_file(self, tmp_path, monkeypatch):
        config_file = tmp_path / "config.json"
        config_file.write_text(json.dumps({
            "server_url": "https://file.malva.bio",
            "api_token": "file-tok",
        }))
        monkeypatch.setenv("MALVA_API_URL", "https://env.malva.bio")

        config = Config.load(str(config_file))
        assert config.server_url == "https://env.malva.bio"
        assert config.api_token == "file-tok"

    def test_load_missing_file(self, tmp_path):
        config = Config.load(str(tmp_path / "missing.json"))
        assert config.server_url is None
        assert config.api_token is None
        assert config.verify_ssl is True


class TestConfigSave:
    def test_save_creates_file(self, tmp_path):
        config_file = tmp_path / "subdir" / "config.json"
        config = Config(str(config_file))
        config.server_url = "https://saved.malva.bio"
        config.api_token = "saved-tok"
        config.save()

        assert config_file.exists()
        data = json.loads(config_file.read_text())
        assert data["server_url"] == "https://saved.malva.bio"
        assert data["api_token"] == "saved-tok"

    def test_save_and_load_roundtrip(self, tmp_path):
        config_file = tmp_path / "config.json"

        c1 = Config(str(config_file))
        c1.server_url = "https://round.malva.bio"
        c1.api_token = "round-tok"
        c1.verify_ssl = False
        c1.save()

        c2 = Config.load(str(config_file))
        assert c2.server_url == c1.server_url
        assert c2.api_token == c1.api_token
        assert c2.verify_ssl == c1.verify_ssl


class TestConfigSetValue:
    def test_aliases(self):
        config = Config()
        config.set_value("server", "https://alias.malva.bio")
        assert config.server_url == "https://alias.malva.bio"

        config.set_value("url", "https://url.malva.bio")
        assert config.server_url == "https://url.malva.bio"

        config.set_value("token", "tok-alias")
        assert config.api_token == "tok-alias"

        config.set_value("ssl", "false")
        assert config.verify_ssl is False

    def test_invalid_key_raises(self):
        config = Config()
        with pytest.raises(ValueError, match="Invalid configuration key"):
            config.set_value("bogus", "value")


class TestConfigValidate:
    def test_missing_url(self):
        config = Config()
        config.api_token = "tok"
        errors = config.validate()
        assert any("server_url" in e for e in errors)

    def test_missing_token(self):
        config = Config()
        config.server_url = "https://valid.malva.bio"
        errors = config.validate()
        assert any("api_token" in e for e in errors)

    def test_valid(self):
        config = Config()
        config.server_url = "https://valid.malva.bio"
        config.api_token = "tok"
        assert config.validate() == []

    def test_invalid_url_scheme(self):
        config = Config()
        config.server_url = "ftp://invalid.malva.bio"
        config.api_token = "tok"
        errors = config.validate()
        assert any("http" in e for e in errors)


class TestConfigShow:
    def test_masks_token(self):
        config = Config()
        config.server_url = "https://show.malva.bio"
        config.api_token = "secret-token"
        shown = config.show()
        assert shown["api_token"] == "***"
        assert shown["server_url"] == "https://show.malva.bio"

    def test_none_token_shown_as_none(self):
        config = Config()
        shown = config.show()
        assert shown["api_token"] is None


class TestConfigReset:
    def test_reset_to_defaults(self):
        config = Config()
        config.server_url = "https://reset.malva.bio"
        config.api_token = "tok"
        config.verify_ssl = False

        config.reset_to_defaults()
        assert config.server_url is None
        assert config.api_token is None
        assert config.verify_ssl is True
